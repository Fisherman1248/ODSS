%% ODSS over a configurable N-path wideband delay-scale channel
% Processing chain:
%
%   x[k,l] --T_iMF--> X[n,m] --ODSS modulator--> s(t)
%          --N-path delay-scale channel + AWGN--> r(t)
%          --ODSS demodulator--> Yhat[n,m]
%          --diagonal one-tap MMSE--> Zhat[n,m]
%          --T_iMF^{-1}--> xhat[k,l]
%
% Channel model:
%   r_s(t) = sum_{p=1}^{P} h_p*sqrt(alpha_p)*
%            s_tx(alpha_p*(t-tau_p))
%
% The channel parameters tau_p, alpha_p and h_p are generated
% automatically from Path_num, tau_max and alpha_max. No path-specific
% delay or scale values are hard-coded.
%
% This script follows the delay-scale-domain model:
%   Yhat = D*X + W
% where D is diagonal and
%   D(i,i) = H_{n_i,m_i}[n_i,m_i].
%
% IMPORTANT:
% The diagonal relation relies on the pulse/lattice assumptions in the
% ODSS paper. The script therefore reports the diagonal-model error:
%   ||Y_signal - D*X|| / ||Y_signal||.
% A large value means the actual sampled implementation has significant
% off-diagonal interference, even though the decoder is intentionally
% using the paper's diagonal approximation.
%
% MATLAB R2016b or newer is required for local functions in scripts.

clear; clc; close all;

%% 1. Reproducibility
dataSeed    = 7;
channelSeed = 21;
noiseSeed   = 31;

%% 2. ODSS parameters
q         = 2.0;
Nscale    = 7;
B         = 1280;                 % total baseband bandwidth, Hz
W         = B*(q-1)/(q^Nscale-1); % basic chirplet bandwidth, Hz
T         = 1.9;                  % basic pulse/block duration, s
Fs        = 10240;                % waveform sampling rate, Hz
Ts        = 1/Fs;
SNRdB     = 30;
pulseType = 'phydyas';            % 'rect' or 'phydyas'

n_set   = 0:Nscale-1;
M_scale = floor(q.^n_set);
M_tot   = sum(M_scale);

%% 3. Configurable random N-path channel
% Change Path_num to any positive integer.
Path_num = 2;

% Path delays are generated uniformly in [0,tau_max].
tau_max = 0.010;                  % maximum delay, seconds

% Path scale factors are generated uniformly in
% [1/alpha_max, alpha_max].
alpha_max = 1.001;                % alpha_max >= 1

% Optional deterministic direct path:
% true  -> path 1 has tau=0 and alpha=1
% false -> all paths are randomly generated
includeDirectPath = true;

% Normalize sum_p |h_p|^2 to one.
normalizePathEnergy = true;

%% 4. Basic parameter validation
validateattributes(Path_num, {'numeric'}, ...
    {'scalar','integer','positive'}, mfilename, 'Path_num');
validateattributes(tau_max, {'numeric'}, ...
    {'scalar','real','nonnegative','finite'}, mfilename, 'tau_max');
validateattributes(alpha_max, {'numeric'}, ...
    {'scalar','real','>=',1,'finite'}, mfilename, 'alpha_max');

if q < alpha_max^2
    warning(['The paper condition q >= alpha_max^2 is not satisfied: ' ...
             'q = %.6g, alpha_max^2 = %.6g.'], q, alpha_max^2);
end

if tau_max > 0
    W_lim_1 = 1 / ((1 + alpha_max)*tau_max);
    W_lim_2 = 1 / ((1 + alpha_max^(2*Nscale-3))*tau_max);
    W_limit = min(W_lim_1, W_lim_2);

    if W > W_limit
        warning(['The selected W may violate the paper''s no-ICI ' ...
                 'condition: W = %.6g Hz, approximate limit = %.6g Hz.'], ...
                 W, W_limit);
    end
else
    W_limit = inf;
end

fprintf('========== ODSS configurable N-path diagonal-MMSE test ==========\n');
fprintf('q = %.3f, Nscale = %d, M_tot = %d\n', q, Nscale, M_tot);
fprintf('B = %.3f Hz, W = %.6f Hz, T = %.6f s, Fs = %.1f Hz\n', ...
    B, W, T, Fs);
fprintf('Path_num = %d, tau_max = %.6f s, alpha_max = %.7f\n', ...
    Path_num, tau_max, alpha_max);
if isfinite(W_limit)
    fprintf('approximate W limit from channel support = %.6f Hz\n\n', ...
        W_limit);
else
    fprintf('\n');
end

%% 5. Scale-major vector index map
% Ordering:
%   (0,0),
%   (1,0),(1,1),
%   (2,0),...,(2,M(2)-1), ...
pair_n = zeros(M_tot,1);
pair_m = zeros(M_tot,1);

idx = 0;
for n = 0:Nscale-1
    for m = 0:M_scale(n+1)-1
        idx = idx + 1;
        pair_n(idx) = n;
        pair_m(idx) = m;
    end
end

%% 6. Discrete inverse Mellin-Fourier transform
% X[n,m] = q^(-n/2)/N *
%          sum_k {1/M(k) * sum_l x[k,l] *
%          exp(j*2*pi*(m*l/M(k)-n*k/N))}
T_iMF = complex(zeros(M_tot, M_tot));

for outIdx = 1:M_tot
    n = pair_n(outIdx);
    m = pair_m(outIdx);

    for inIdx = 1:M_tot
        k  = pair_n(inIdx);
        l  = pair_m(inIdx);
        Mk = M_scale(k+1);

        T_iMF(outIdx,inIdx) = ...
            q^(-n/2) / (Nscale*Mk) * ...
            exp(1j*2*pi*(m*l/Mk - n*k/Nscale));
    end
end

fprintf('Transform diagnostics\n');
fprintf('rank(T_iMF) = %d / %d\n', rank(T_iMF), M_tot);
fprintf('cond(T_iMF) = %.3e\n', cond(T_iMF));
fprintf('inverse check = %.3e\n\n', ...
    norm((T_iMF\T_iMF)-eye(M_tot),'fro')/sqrt(M_tot));

%% 7. Build sampled ODSS transmit basis
Ns = round(T*Fs);
t  = (0:Ns-1).' / Fs;

G_tx = complex(zeros(Ns, M_tot));

for col = 1:M_tot
    n = pair_n(col);
    m = pair_m(col);

    tLocal = q^n .* (t - m/(q^n*W));

    G_tx(:,col) = q^(n/2) .* ...
        odssTransmitPulse(tLocal, q, W, T, pulseType);
end

%% 8. Canonical-dual receive basis for the AWGN smoke-test stage
Gram_tx = Ts * (G_tx' * G_tx);

if rcond(Gram_tx) < 1e-12
    warning('The transmit Gram matrix is close to singular.');
end

G_rx = G_tx / Gram_tx;

biorthErr = norm( ...
    Ts*(G_rx'*G_tx)-eye(M_tot), 'fro') / sqrt(M_tot);

fprintf('Waveform diagnostics\n');
fprintf('samples = %d\n', Ns);
fprintf('rank(G_tx) = %d / %d\n', rank(G_tx), M_tot);
fprintf('cond(Ts*G_tx^H*G_tx) = %.3e\n', cond(Gram_tx));
fprintf('biorthogonality error = %.3e\n\n', biorthErr);

%% 9. Generate unit-power Gray QPSK symbols
rng(dataSeed, 'twister');

bits = randi([0 1], 2*M_tot, 1);

x_vector = ((1 - 2*bits(1:2:end)) + ...
         1j*(1 - 2*bits(2:2:end))) / sqrt(2);

%% 10. ODSS transform and waveform modulation
X_vector = T_iMF * x_vector;
s_tx     = G_tx * X_vector;

%% 11. Noiseless identity-channel closed-loop check
Y0_vector = Ts * (G_rx' * s_tx);
x0_soft   = T_iMF \ Y0_vector;

fprintf('Identity-channel closed-loop\n');
fprintf('||Y0-X||/||X|| = %.3e\n', ...
    norm(Y0_vector-X_vector)/max(norm(X_vector),eps));
fprintf('||x0-x||/||x|| = %.3e\n\n', ...
    norm(x0_soft-x_vector)/max(norm(x_vector),eps));

%% 12. Generate random N-path channel parameters
rng(channelSeed, 'twister');

tau_p = tau_max * rand(Path_num,1);

alpha_min = 1/alpha_max;
alpha_p = alpha_min + (alpha_max-alpha_min)*rand(Path_num,1);

h_p = (randn(Path_num,1) + 1j*randn(Path_num,1))/sqrt(2);

if includeDirectPath
    tau_p(1)   = 0;
    alpha_p(1) = 1;
end

if normalizePathEnergy
    hNorm = norm(h_p);
    if hNorm == 0
        error('Generated path-gain vector has zero norm.');
    end
    h_p = h_p/hNorm;
end

% Sort paths by increasing delay to make printed output easier to read.
[tau_p, sortIdx] = sort(tau_p);
alpha_p = alpha_p(sortIdx);
h_p     = h_p(sortIdx);

fprintf('Generated channel paths\n');
for p = 1:Path_num
    fprintf(['path %3d: tau = %9.6f s, alpha = %.7f, ' ...
             'h = %+.6f %+.6fj, |h| = %.6f\n'], ...
        p, tau_p(p), alpha_p(p), ...
        real(h_p(p)), imag(h_p(p)), abs(h_p(p)));
end
fprintf('sum |h_p|^2 = %.6f\n\n', sum(abs(h_p).^2));

%% 13. Construct the received time axis
% A scaled pulse s(alpha*(t-tau)) has duration T/alpha and ends at
% tau + T/alpha.
tEnd = max([T; tau_p + T./alpha_p]);
Nrx  = ceil(tEnd*Fs) + 1;
t_rx = (0:Nrx-1).' / Fs;

%% 14. Pass the complete ODSS waveform through the N-path channel
% r_s(t) = sum_p h_p*sqrt(alpha_p)*s_tx(alpha_p*(t-tau_p))
r_s = complex(zeros(Nrx,1));

for p = 1:Path_num
    tQuery = alpha_p(p) .* (t_rx - tau_p(p));

    s_path = interp1( ...
        t, s_tx, tQuery, 'linear', 0);

    r_s = r_s + ...
        h_p(p)*sqrt(alpha_p(p))*s_path;
end

%% 15. Add time-domain complex AWGN
signalPower = mean(abs(r_s).^2);
noisePower  = signalPower / 10^(SNRdB/10);

rng(noiseSeed, 'twister');
w_time = sqrt(noisePower/2) .* ...
    (randn(Nrx,1) + 1j*randn(Nrx,1));

r = r_s + w_time;

measuredSNRdB = 10*log10( ...
    mean(abs(r_s).^2)/mean(abs(w_time).^2));

fprintf('Received waveform\n');
fprintf('receiver duration = %.6f s\n', Nrx/Fs);
fprintf('received samples  = %d\n', Nrx);
fprintf('signal power      = %.6e\n', signalPower);
fprintf('noise power       = %.6e\n', noisePower);
fprintf('target SNR        = %.2f dB\n', SNRdB);
fprintf('measured SNR      = %.2f dB\n\n', measuredSNRdB);

%% 16. Extend the receive basis to the received time axis
% The receive waveforms are zero outside their original sampled support.
G_rx_ext = complex(zeros(Nrx, M_tot));
G_rx_ext(1:Ns,:) = G_rx;

%% 17. ODSS demodulator outputs
Yhat_vector = Ts * (G_rx_ext' * r);
Ysig_vector = Ts * (G_rx_ext' * r_s);
W_vector    = Ts * (G_rx_ext' * w_time);

demodLinearityErr = norm( ...
    Yhat_vector-(Ysig_vector+W_vector)) / ...
    max(norm(Yhat_vector),eps);

fprintf('ODSS demodulator\n');
fprintf('||Yhat-(Ysignal+W)||/||Yhat|| = %.3e\n\n', ...
    demodLinearityErr);

%% 18. Compute only the diagonal delay-scale channel D
% No M_tot-by-M_tot effective channel matrix is constructed.
%
% For coefficient i corresponding to lattice point (n_i,m_i):
%   d_i = H_{n_i,m_i}[n_i,m_i]
%
% Each path's contribution is obtained by passing every transmit basis
% waveform through that path and taking only the same-index receive inner
% products.

d_vector = complex(zeros(M_tot,1));

for p = 1:Path_num
    tQuery = alpha_p(p) .* (t_rx - tau_p(p));

    % Channel-warp all transmit basis functions for path p.
    % G_path has size Nrx-by-M_tot.
    G_path = interp1( ...
        t, G_tx, tQuery, 'linear', 0);

    % Same-index inner products only:
    % d_i += h_p*sqrt(alpha_p)*<g_rx_i, channel(g_tx_i)>.
    sameIndexInnerProducts = ...
        Ts * sum(conj(G_rx_ext).*G_path, 1).';

    d_vector = d_vector + ...
        h_p(p)*sqrt(alpha_p(p))*sameIndexInnerProducts;
end

% D is diagonal by construction in the paper's delay-scale decoder.
D = diag(d_vector); %#ok<NASGU>

%% 19. Check the accuracy of the diagonal model
Ydiag_vector = d_vector .* X_vector;
ICI_ISI_vector = Ysig_vector - Ydiag_vector;

diagModelErr = norm(ICI_ISI_vector) / ...
               max(norm(Ysig_vector),eps);

desiredPowerDS = mean(abs(Ydiag_vector).^2);
interfPowerDS  = mean(abs(ICI_ISI_vector).^2);

if interfPowerDS > 0
    diagModelSIRdB = 10*log10(desiredPowerDS/interfPowerDS);
else
    diagModelSIRdB = inf;
end

fprintf('Diagonal delay-scale model\n');
fprintf('min |d_i| = %.6e\n', min(abs(d_vector)));
fprintf('max |d_i| = %.6e\n', max(abs(d_vector)));
fprintf('mean |d_i| = %.6e\n', mean(abs(d_vector)));
fprintf('||Ysignal-DX||/||Ysignal|| = %.3e\n', diagModelErr);
fprintf('diagonal-model desired/interference ratio = %.2f dB\n\n', ...
    diagModelSIRdB);

%% 20. Delay-scale-domain one-tap MMSE equalizer
% The paper writes:
%   Zhat = D^H*(D*D^H + sigma_W^2*I)^(-1)*Yhat.
%
% Since D is diagonal, implement it element-by-element.
%
% The current discrete iMF implementation does not produce unit-power X
% under ordinary vector stacking, so use a power-normalized loading:
%   lambda = average delay-scale noise variance / average X power.
%
% Compute the average output-noise variance analytically from the known
% time-domain AWGN variance instead of estimating it from this one random
% realization.

noiseVarDS_each = noisePower * Ts^2 .* ...
    sum(abs(G_rx_ext).^2, 1).';

sigmaW2_delayScale = mean(noiseVarDS_each);
Es_delayScale      = mean(abs(X_vector).^2);
mmseLoading        = sigmaW2_delayScale / max(Es_delayScale,eps);

Zhat_vector = conj(d_vector) ./ ...
    (abs(d_vector).^2 + mmseLoading) .* Yhat_vector;

fprintf('MMSE parameters\n');
fprintf('delay-scale signal power = %.6e\n', Es_delayScale);
fprintf('average delay-scale noise variance = %.6e\n', ...
    sigmaW2_delayScale);
fprintf('power-normalized MMSE loading = %.6e\n\n', mmseLoading);

%% 21. Inverse iMF and QPSK slicing
x_soft = T_iMF \ Zhat_vector;

bits_hat = zeros(size(bits));
bits_hat(1:2:end) = real(x_soft) < 0;
bits_hat(2:2:end) = imag(x_soft) < 0;

x_hat = ((1 - 2*bits_hat(1:2:end)) + ...
      1j*(1 - 2*bits_hat(2:2:end))) / sqrt(2);

numBitErrors = nnz(bits_hat ~= bits);
BER          = numBitErrors/numel(bits);

fprintf('Recovered symbols\n');
fprintf('||Zhat-X||/||X|| = %.3e\n', ...
    norm(Zhat_vector-X_vector)/max(norm(X_vector),eps));
fprintf('||xsoft-x||/||x|| = %.3e\n', ...
    norm(x_soft-x_vector)/max(norm(x_vector),eps));
fprintf('bit errors = %d / %d\n', numBitErrors, numel(bits));
fprintf('BER = %.6e\n', BER);

%% 22. Plots
figure;
plot(real(x_soft), imag(x_soft), 'o');
hold on;
plot(real(x_vector), imag(x_vector), 'kx', 'LineWidth', 1.5);
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title(sprintf( ...
    'ODSS, %d-path diagonal MMSE, SNR = %.1f dB', ...
    Path_num, SNRdB));
legend('Before slicing', 'Transmitted QPSK', ...
    'Location', 'best');

figure;
stem(1:M_tot, abs(d_vector), 'filled');
grid on;
xlabel('Stacked delay-scale coefficient index');
ylabel('|d_i|');
title('Diagonal delay-scale channel magnitudes');


%% Local function: basic chirplet and optional PHYDYAS window
function g_tx = odssTransmitPulse(tLocal, q, W, T, pulseType)
% Basic linear chirplet with physical bandwidth W and support [0,T).

    g_tx = complex(zeros(size(tLocal)));

    valid = (tLocal >= 0) & (tLocal < T);
    if ~any(valid)
        return;
    end

    tv = tLocal(valid);

    freqScale = W/(sqrt(q)-1/sqrt(q));
    f1 = freqScale/sqrt(q);
    f2 = freqScale*sqrt(q);
    kappa = (f2-f1)/T;

    g0 = exp(1j*2*pi*(f1*tv + 0.5*kappa*tv.^2));

    switch lower(pulseType)
        case 'rect'
            window = ones(size(tv));

        case 'phydyas'
            K = 3;
            A = [0.91143783, 0.41143783];

            window = ones(size(tv));
            for kk = 1:K-1
                window = window + ...
                    2*(-1)^kk*A(kk).* ...
                    cos(2*pi*kk*tv/(K*T));
            end

        otherwise
            error('Unknown pulseType: %s', pulseType);
    end

    g_tx(valid) = window .* g0;
end
