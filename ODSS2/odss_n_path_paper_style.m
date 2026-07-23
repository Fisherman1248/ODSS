%% ODSS paper-style N-path simulation
% Paper:
% A. K. P. and C. R. Murthy,
% "Orthogonal Delay Scale Space Modulation: A New Technique for
% Wideband Time-Varying Channels," IEEE TSP, 2022.
%
% This script follows the paper more closely than the earlier smoke test:
%   1) PHYDYAS-windowed chirplet, Eqs. (71)-(73)
%   2) matched-filter receiver (G_rx = G_tx), not a canonical dual
%   3) BPSK, as in the BER simulations
%   4) configurable P-path channel, h_p ~ CN(0,1)
%   5) tau_p ~ U(0,tau_max)
%   6) alpha_p ~ U(1/alpha_max,alpha_max)
%   7) 8x channel oversampling
%   8) delays rounded to the 8x time grid
%   9) rational resampling approximation of alpha_p
%  10) downsampling back to the original waveform rate
%  11) delay-scale one-tap MMSE using only diag(D)
%
% IMPORTANT IMPLEMENTATION NOTE
% The paper prints normalized chirplet endpoints f1=1/sqrt(q) and
% f2=sqrt(q), but does not explicitly state the dimensional frequency
% scale used to generate Fig. 8. Here they are mapped to physical
% frequencies so that f2-f1 = W:
%
%       f1 = W/(q-1),  f2 = qW/(q-1).
%
% For q=2 this gives f1=W and f2=2W.
%
% The script offers two lattice/window interpretations:
%
%   cfg.latticeMode = 'equation73'
%       delay = m/(q^n W), exactly as printed in Eq. (73).
%
%   cfg.latticeMode = 'figure14'
%       delay = m*T/M(n), which visibly tiles M(n) pulses across the
%       T-second block as in the supplementary waveform figure.
%
%   cfg.phydyasMode = 'literal'
%       uses cos(2*pi*k*t/(K*T)), exactly as printed in Eq. (72).
%
%   cfg.phydyasMode = 'symmetric'
%       uses a centered finite-duration PHYDYAS prototype by interpreting
%       the printed K*T denominator as the complete prototype duration.
%
% Start with the literal equation settings. The alternatives are included
% because the paper does not provide all dimensional/discrete-time details
% needed for a bit-for-bit reproduction of every figure.
%
% Requires Signal Processing Toolbox for resample().
% MATLAB R2016b or newer is required for local functions in a script.

clear;
clc;
close all;

%% 1. User configuration
cfg.dataSeed       = 7;
cfg.channelSeed    = 21;
cfg.noiseSeed      = 31;

% ODSS parameters used in the paper example
cfg.q              = 2;
cfg.Nscale         = 7;          % n = 0,...,6
cfg.B              = 1280;       % downconverted baseband width, Hz
cfg.T              = 1.9;        % ODSS symbol-block duration, seconds
cfg.Fs             = 1280;       % receiver/base waveform rate, Hz
cfg.channelOS      = 8;          % paper: oversample channel by 8
cfg.SNRdB          = 30;

% Configurable N-path channel
cfg.Path_num       = 2;          % change to any positive integer
cfg.tau_max        = 10e-3;
cfg.alpha_max      = 1.001;

% Paper-style channel generation
cfg.forceDirectPath       = false;
cfg.normalizePathEnergy   = false;
cfg.alphaApproxTolerance  = 1e-5;

% Pulse/lattice interpretation
cfg.pulseType       = 'phydyas';   % 'phydyas' or 'rect'
cfg.latticeMode     = 'equation73';% 'equation73' or 'figure14'
cfg.phydyasMode     = 'literal';   % 'literal' or 'symmetric'

% Receiver/equalizer
cfg.mmseMode        = 'paper';     % 'paper' or 'power_normalized'
cfg.computeFullDiagnostic = true;  % full matrix is diagnostic only

% Plot controls
cfg.makePlots       = true;

%% 2. Derived parameters and validation
if exist('resample','file') ~= 2
    error(['This paper-style channel implementation requires resample(), ' ...
           'which is provided by Signal Processing Toolbox.']);
end

validateattributes(cfg.Path_num, {'numeric'}, ...
    {'scalar','integer','positive'}, mfilename, 'cfg.Path_num');
validateattributes(cfg.alpha_max, {'numeric'}, ...
    {'scalar','real','finite','>=',1}, mfilename, 'cfg.alpha_max');
validateattributes(cfg.tau_max, {'numeric'}, ...
    {'scalar','real','finite','nonnegative'}, mfilename, 'cfg.tau_max');

q        = cfg.q;
Nscale   = cfg.Nscale;
B        = cfg.B;
T        = cfg.T;
Fs       = cfg.Fs;
Ts       = 1/Fs;
FsHi     = cfg.channelOS*Fs;

% Eq. (77)
W = B*(q-1)/(q^Nscale-1);

n_set   = 0:Nscale-1;
M_scale = ceil(q.^n_set);
M_tot   = sum(M_scale);

fprintf('=============== ODSS paper-style N-path test ===============\n');
fprintf('q = %.3f, Nscale = %d, M_tot = %d\n', q, Nscale, M_tot);
fprintf('B = %.3f Hz, W = %.9f Hz, T = %.6f s\n', B, W, T);
fprintf('Fs = %.1f Hz, channel Fs = %.1f Hz (%dx)\n', ...
    Fs, FsHi, cfg.channelOS);
fprintf('pulse = %s, lattice = %s, PHYDYAS mode = %s\n', ...
    cfg.pulseType, cfg.latticeMode, cfg.phydyasMode);
fprintf('Path_num = %d, tau_max = %.6f s, alpha_max = %.7f\n\n', ...
    cfg.Path_num, cfg.tau_max, cfg.alpha_max);

if q < cfg.alpha_max^2
    warning('Paper condition q >= alpha_max^2 is not satisfied.');
end

if cfg.tau_max > 0
    W_limit_1 = 1/((1+cfg.alpha_max)*cfg.tau_max);
    W_limit_2 = 1/((1+cfg.alpha_max^(2*Nscale-3))*cfg.tau_max);
    W_limit   = min(W_limit_1,W_limit_2);

    fprintf('Approximate paper W limit = %.6f Hz\n', W_limit);
    fprintf('Selected W                = %.6f Hz\n\n', W);

    if W > W_limit
        warning('Selected W exceeds the approximate no-ICI support limit.');
    end
end

%% 3. Scale-major index map
pair_n = zeros(M_tot,1);
pair_m = zeros(M_tot,1);

index = 0;
for n = 0:Nscale-1
    for m = 0:M_scale(n+1)-1
        index = index+1;
        pair_n(index) = n;
        pair_m(index) = m;
    end
end

%% 4. Discrete inverse Mellin-Fourier transform, Eq. (38)
T_iMF = complex(zeros(M_tot,M_tot));

for outIndex = 1:M_tot
    n = pair_n(outIndex);
    m = pair_m(outIndex);

    for inIndex = 1:M_tot
        k  = pair_n(inIndex);
        l  = pair_m(inIndex);
        Mk = M_scale(k+1);

        T_iMF(outIndex,inIndex) = ...
            q^(-n/2)/(Nscale*Mk) * ...
            exp(1j*2*pi*(m*l/Mk-n*k/Nscale));
    end
end

fprintf('Transform diagnostics\n');
fprintf('rank(T_iMF) = %d / %d\n', rank(T_iMF), M_tot);
fprintf('cond(T_iMF) = %.3e\n', cond(T_iMF));
fprintf('inverse check = %.3e\n', ...
    norm((T_iMF\T_iMF)-eye(M_tot),'fro')/sqrt(M_tot));
fprintf('Euclidean energy-preservation error = %.3e\n\n', ...
    norm(T_iMF'*T_iMF-eye(M_tot),'fro')/sqrt(M_tot));

%% 5. Build the sampled PHYDYAS ODSS transmit basis
Ns = round(T*Fs);
t  = (0:Ns-1).'/Fs;

% Physical mapping of the normalized endpoints in Eq. (71).
f1 = W/(q-1);
f2 = q*W/(q-1);
chirpRate = (f2-f1)/T;

% Normalize the prototype to unit sampled continuous-time energy.
gPrototype = odssPrototypePulse( ...
    t, T, f1, chirpRate, cfg.pulseType, cfg.phydyasMode);

prototypeEnergy = Ts*sum(abs(gPrototype).^2);
prototypeNorm   = sqrt(prototypeEnergy);

if prototypeNorm == 0
    error('The generated prototype pulse has zero energy.');
end

G_tx = complex(zeros(Ns,M_tot));

for column = 1:M_tot
    n  = pair_n(column);
    m  = pair_m(column);
    Mn = M_scale(n+1);

    switch lower(cfg.latticeMode)
        case 'equation73'
            delay = m/(q^n*W);

        case 'figure14'
            delay = m*T/Mn;

        otherwise
            error('Unsupported cfg.latticeMode: %s',cfg.latticeMode);
    end

    localTime = q^n*(t-delay);

    G_tx(:,column) = q^(n/2) * ...
        odssPrototypePulse( ...
        localTime,T,f1,chirpRate,cfg.pulseType,cfg.phydyasMode) ...
        /prototypeNorm;
end

%% 6. Paper-style receive pulse: matched filtering
% No canonical-dual matrix inversion is used here.
G_rx = G_tx;

Gram = Ts*(G_rx'*G_tx);
columnEnergy = real(diag(Gram));
GramNormalized = Gram ./ sqrt(columnEnergy*columnEnergy.');
GramOffDiagonal = GramNormalized-eye(M_tot);

maxRawCorrelation = max(abs(GramOffDiagonal(:)));

fprintf('Matched-filter waveform diagnostics\n');
fprintf('samples = %d\n', Ns);
fprintf('rank(G_tx) = %d / %d\n', rank(G_tx), M_tot);
fprintf('min/max column energy = %.6e / %.6e\n', ...
    min(columnEnergy),max(columnEnergy));
fprintf('cond(Ts*G_tx^H*G_tx) = %.3e\n',cond(Gram));
fprintf('normalized Gram error = %.3e\n', ...
    norm(GramOffDiagonal,'fro')/sqrt(M_tot));
fprintf('maximum pairwise correlation = %.6e (%.2f dB)\n\n', ...
    maxRawCorrelation,20*log10(max(maxRawCorrelation,realmin)));

%% 7. BPSK data, as in the paper BER simulations
rng(cfg.dataSeed,'twister');

bits = randi([0 1],M_tot,1);
x_vector = 2*bits-1;

X_vector = T_iMF*x_vector;
s_tx     = G_tx*X_vector;

fprintf('Data and waveform power\n');
fprintf('mean |x|^2 = %.6e\n',mean(abs(x_vector).^2));
fprintf('mean |X|^2 = %.6e\n',mean(abs(X_vector).^2));
fprintf('mean |s|^2 = %.6e\n\n',mean(abs(s_tx).^2));

%% 8. Identity-channel matched-filter check
Y_identity = Ts*(G_rx'*s_tx);
d_identity = diag(Gram);

% Diagonal equalization of the identity-channel matched-filter output.
Z_identity = Y_identity./d_identity;
x_identity = T_iMF\Z_identity;

fprintf('Identity-channel matched-filter check\n');
fprintf('||Y-GX||/||Y|| = %.3e\n', ...
    norm(Y_identity-Gram*X_vector)/max(norm(Y_identity),eps));
fprintf('diagonal-model error = %.3e\n', ...
    norm(Y_identity-d_identity.*X_vector)/max(norm(Y_identity),eps));
fprintf('||x_identity-x||/||x|| = %.3e\n\n', ...
    norm(x_identity-x_vector)/norm(x_vector));

%% 9. Generate the random P-path channel
rng(cfg.channelSeed,'twister');

tau_p = cfg.tau_max*rand(cfg.Path_num,1);

alphaLower = 1/cfg.alpha_max;
alpha_p = alphaLower + ...
    (cfg.alpha_max-alphaLower)*rand(cfg.Path_num,1);

h_p = (randn(cfg.Path_num,1)+1j*randn(cfg.Path_num,1))/sqrt(2);

if cfg.forceDirectPath
    tau_p(1)   = 0;
    alpha_p(1) = 1;
end

if cfg.normalizePathEnergy
    h_p = h_p/norm(h_p);
end

% Sort only for readable output.
[tau_p,sortOrder] = sort(tau_p);
alpha_p = alpha_p(sortOrder);
h_p     = h_p(sortOrder);

fprintf('Generated channel paths\n');
for p = 1:cfg.Path_num
    fprintf(['path %3d: tau = %.9f s, alpha = %.9f, ' ...
             'h = %+.6f %+.6fj, |h| = %.6f\n'], ...
        p,tau_p(p),alpha_p(p),real(h_p(p)),imag(h_p(p)),abs(h_p(p)));
end
fprintf('sum |h_p|^2 = %.6f\n\n',sum(abs(h_p).^2));

%% 10. Paper-style 8x oversampled delay-scale channel
% Apply the same channel to every transmit basis column. This gives an
% exact linear sampled channel representation:
%
%       r_signal = G_channel * X.
%
% The full effective delay-scale matrix is used only for diagnostics.
[G_channel,channelInfo] = applyPaperDelayScaleChannel( ...
    G_tx,Fs,cfg.channelOS,h_p,tau_p,alpha_p, ...
    cfg.alphaApproxTolerance);

Nrx = size(G_channel,1);
t_rx = (0:Nrx-1).'/Fs;

r_signal = G_channel*X_vector;

fprintf('Rational channel approximations\n');
for p = 1:cfg.Path_num
    fprintf(['path %3d: tau requested/used = %.9f / %.9f s, ' ...
             'alpha requested/used = %.9f / %.9f, error = %.3e\n'], ...
        p,tau_p(p),channelInfo.tauUsed(p), ...
        alpha_p(p),channelInfo.alphaUsed(p), ...
        abs(alpha_p(p)-channelInfo.alphaUsed(p)));
end
fprintf('\n');

%% 11. Add AWGN at the original receiver waveform rate
signalPower = mean(abs(r_signal).^2);
noisePower  = signalPower/10^(cfg.SNRdB/10);

rng(cfg.noiseSeed,'twister');
w_time = sqrt(noisePower/2) * ...
    (randn(Nrx,1)+1j*randn(Nrx,1));

r = r_signal+w_time;

measuredSNRdB = 10*log10( ...
    mean(abs(r_signal).^2)/mean(abs(w_time).^2));

fprintf('Received waveform\n');
fprintf('receiver duration = %.6f s\n',Nrx/Fs);
fprintf('received samples = %d\n',Nrx);
fprintf('signal power = %.6e\n',signalPower);
fprintf('noise power = %.6e\n',noisePower);
fprintf('target/measured SNR = %.2f / %.2f dB\n\n', ...
    cfg.SNRdB,measuredSNRdB);

%% 12. Matched-filter ODSS demodulator
G_rx_ext = complex(zeros(Nrx,M_tot));
copyLength = min(Ns,Nrx);
G_rx_ext(1:copyLength,:) = G_rx(1:copyLength,:);

Yhat_vector = Ts*(G_rx_ext'*r);
Ysig_vector = Ts*(G_rx_ext'*r_signal);
W_vector    = Ts*(G_rx_ext'*w_time);

fprintf('Demodulator linearity\n');
fprintf('||Yhat-(Ysignal+W)||/||Yhat|| = %.3e\n\n', ...
    norm(Yhat_vector-(Ysig_vector+W_vector))/ ...
    max(norm(Yhat_vector),eps));

%% 13. Diagonal D used by the paper one-tap receiver
% Same-index matched-filter outputs only:
d_vector = Ts*sum(conj(G_rx_ext).*G_channel,1).';

Ydiag_vector = d_vector.*X_vector;
ICI_vector   = Ysig_vector-Ydiag_vector;

diagModelError = norm(ICI_vector)/max(norm(Ysig_vector),eps);
desiredPower   = mean(abs(Ydiag_vector).^2);
iciPower       = mean(abs(ICI_vector).^2);

if iciPower > 0
    diagonalSIRdB = 10*log10(desiredPower/iciPower);
else
    diagonalSIRdB = inf;
end

fprintf('Delay-scale diagonal model\n');
fprintf('min/mean/max |d_i| = %.6e / %.6e / %.6e\n', ...
    min(abs(d_vector)),mean(abs(d_vector)),max(abs(d_vector)));
fprintf('||Ysignal-DX||/||Ysignal|| = %.3e\n',diagModelError);
fprintf('desired/ICI ratio = %.2f dB\n',diagonalSIRdB);

if cfg.computeFullDiagnostic
    H_DS = Ts*(G_rx_ext'*G_channel);
    D_diagnostic = diag(diag(H_DS));

    offDiagonalRatio = norm(H_DS-D_diagnostic,'fro') / ...
        max(norm(H_DS,'fro'),eps);

    rowICI = sum(abs(H_DS-D_diagnostic).^2,2);
    rowDesired = abs(diag(H_DS)).^2;
    maxICIdB = max(10*log10(max(rowICI,realmin)./ ...
        max(rowDesired,realmin)));

    fprintf('full diagnostic off-diagonal ratio = %.3e\n', ...
        offDiagonalRatio);
    fprintf('maximum row ICI/desired ratio = %.2f dB\n',maxICIdB);
end
fprintf('\n');

%% 14. One-tap MMSE, Eq. (69)
% Exact matched-filter output-noise variance for each coefficient:
noiseVarianceEach = noisePower*Ts^2 * ...
    sum(abs(G_rx_ext).^2,1).';

sigmaW2 = mean(noiseVarianceEach);

switch lower(cfg.mmseMode)
    case 'paper'
        % Literal scalar loading in Eq. (69).
        mmseLoading = sigmaW2;

    case 'power_normalized'
        % Useful diagnostic when the implemented vector form of Eq. (38)
        % is not Euclidean-unitary.
        EsX = mean(abs(X_vector).^2);
        mmseLoading = sigmaW2/max(EsX,eps);

    otherwise
        error('Unsupported cfg.mmseMode: %s',cfg.mmseMode);
end

Zhat_vector = conj(d_vector) ./ ...
    (abs(d_vector).^2+mmseLoading) .* Yhat_vector;

fprintf('MMSE parameters\n');
fprintf('mean delay-scale noise variance = %.6e\n',sigmaW2);
fprintf('loading used = %.6e (%s mode)\n\n', ...
    mmseLoading,cfg.mmseMode);

%% 15. Inverse ODSS transform and BPSK slicing
x_soft = T_iMF\Zhat_vector;
bits_hat = real(x_soft) >= 0;
x_hat = 2*bits_hat-1;

numBitErrors = nnz(bits_hat~=bits);
BER = numBitErrors/numel(bits);

fprintf('Recovered data\n');
fprintf('||Zhat-X||/||X|| = %.3e\n', ...
    norm(Zhat_vector-X_vector)/max(norm(X_vector),eps));
fprintf('||xsoft-x||/||x|| = %.3e\n', ...
    norm(x_soft-x_vector)/norm(x_vector));
fprintf('bit errors = %d / %d\n',numBitErrors,numel(bits));
fprintf('BER = %.6e\n',BER);

%% 16. Plots
if cfg.makePlots
    figure('Name','ODSS matched-filter Gram matrix');
    imagesc(20*log10(abs(GramNormalized)+eps));
    axis image;
    colorbar;
    caxis([-120 0]);
    xlabel('subcarrier index');
    ylabel('subcarrier index');
    title('Normalized transmit-waveform correlations (dB)');

    figure('Name','ODSS delay-scale channel matrix');
    if cfg.computeFullDiagnostic
        Hplot = H_DS ./ max(abs(diag(H_DS)).',eps);
        imagesc(20*log10(abs(Hplot)+eps));
        axis image;
        colorbar;
        caxis([-80 10]);
        xlabel('transmit coefficient index');
        ylabel('receive coefficient index');
        title('Delay-scale effective channel, row-normalized magnitude (dB)');
    else
        stem(1:M_tot,abs(d_vector),'filled');
        grid on;
        xlabel('coefficient index');
        ylabel('|d_i|');
        title('Diagonal delay-scale gains');
    end

    figure('Name','Recovered BPSK symbols');
    plot(real(x_soft),imag(x_soft),'o');
    hold on;
    plot(real(x_vector),imag(x_vector),'kx','LineWidth',1.5);
    grid on;
    xlabel('real');
    ylabel('imaginary');
    title(sprintf('ODSS paper-style receiver, P=%d, SNR=%.1f dB', ...
        cfg.Path_num,cfg.SNRdB));
    legend('soft estimates','transmitted BPSK','Location','best');
end


%% Local function: Eq. (71) chirplet and Eq. (72) PHYDYAS window
function pulse = odssPrototypePulse( ...
    time,T,f1,chirpRate,pulseType,phydyasMode)

    pulse = complex(zeros(size(time)));

    support = (time>=0)&(time<T);
    if ~any(support)
        return;
    end

    localTime = time(support);

    chirplet = exp(1j*2*pi*( ...
        f1*localTime+0.5*chirpRate*localTime.^2));

    switch lower(pulseType)
        case 'rect'
            window = ones(size(localTime));

        case 'phydyas'
            K  = 3;
            Ak = [0.91143783 0.41143783];

            switch lower(phydyasMode)
                case 'literal'
                    % Eq. (72), exactly as printed.
                    denominator = K*T;

                case 'symmetric'
                    % Finite-duration centered PHYDYAS interpretation.
                    denominator = T;

                otherwise
                    error('Unsupported PHYDYAS mode: %s',phydyasMode);
            end

            window = ones(size(localTime));
            for k = 1:K-1
                window = window + ...
                    2*(-1)^k*Ak(k).* ...
                    cos(2*pi*k*localTime/denominator);
            end

        otherwise
            error('Unsupported pulse type: %s',pulseType);
    end

    pulse(support) = window.*chirplet;
end


%% Local function: paper-style oversampled delay-scale channel
function [outputBase,info] = applyPaperDelayScaleChannel( ...
    inputBase,Fs,oversampling,h,tau,alpha,alphaTolerance)
%
% Implements:
%   y(t) = sum_p h_p*sqrt(alpha_p)*x(alpha_p*(t-tau_p))
%
% Procedure follows the simulation description in the paper:
%   - oversample by oversampling
%   - round tau_p to the high-rate grid
%   - rationally approximate alpha_p
%   - resample each path
%   - add the paths
%   - downsample to the original rate

    inputHigh = resample(inputBase,oversampling,1);
    FsHigh = oversampling*Fs;

    numberPaths = numel(h);
    numberColumns = size(inputBase,2);

    outputHigh = complex(zeros(0,numberColumns));

    tauUsed   = zeros(numberPaths,1);
    alphaUsed = zeros(numberPaths,1);

    for p = 1:numberPaths
        delaySamples = round(tau(p)*FsHigh);
        tauUsed(p) = delaySamples/FsHigh;

        % numerator/denominator approximates alpha.
        [alphaNumerator,alphaDenominator] = ...
            rat(alpha(p),alphaTolerance);

        alphaUsed(p) = alphaNumerator/alphaDenominator;

        % resample(x,q,p) gives x(alpha*t) when p/q ~= alpha.
        scaledPath = resample( ...
            inputHigh,alphaDenominator,alphaNumerator);

        pathLength = delaySamples+size(scaledPath,1);

        if size(outputHigh,1)<pathLength
            outputHigh(pathLength,numberColumns) = 0;
        end

        rowRange = delaySamples+(1:size(scaledPath,1));

        outputHigh(rowRange,:) = outputHigh(rowRange,:) + ...
            h(p)*sqrt(alphaUsed(p))*scaledPath;
    end

    outputBase = resample(outputHigh,1,oversampling);

    info.tauUsed   = tauUsed;
    info.alphaUsed = alphaUsed;
end
