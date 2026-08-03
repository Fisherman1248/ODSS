%% ODSS over a pure AWGN channel
% Follows the slide/paper processing chain:
%
%   x[k,l] --Eq. (38)--> X[n,m] --ODSS modulator--> s(t)
%          --AWGN--> r(t) --ODSS demodulator--> Yhat[n,m]
%          --one-tap MMSE--> Zhat[n,m] --T_iMF^{-1}--> xhat[k,l]
%
% Pure AWGN channel:
%   h(tau,alpha) = delta(tau) delta(alpha-1)
%   D = I
%   Yhat = X + W
%
% The sampled receiver below uses the canonical dual basis so that
% Ts*G_rx'*G_tx = I numerically.
%
% MATLAB R2016b or newer is required for local functions in scripts.

clear; clc; close all;

%% 1. Parameters
q = 2.0;
Nscale = 6;
B = 1280;
W = B*(q-1)/(q^Nscale-1);
T = 1.9;
Fs = 10240;
Ts = 1/Fs;
SNRdB = 20;
pulseType = 'phydyas';
QAM_order = 4;

n_set = 0:Nscale-1;
M_scale = floor(q.^n_set);
M_tot = sum(M_scale);

fprintf('================ ODSS pure-AWGN test ================\n');
fprintf('q = %.3f, Nscale = %d, M_tot = %d\n', q, Nscale, M_tot);
fprintf('B = %.3f Hz, W = %.6f Hz, T = %.6f s, Fs = %.1f Hz\n\n', B, W, T, Fs);

%% 2. Valid ODSS index mapping
% pair_n(idx), pair_m(idx) give the mathematical index (n,m)
% associated with the idx-th element of X_vector.
pair_n = zeros(M_tot,1);
pair_m = zeros(M_tot,1);

indexPointer = 1;

for n = 0:Nscale-1
    Mn = M_scale(n+1);

    for m = 0:Mn-1
        pair_n(indexPointer) = n;
        pair_m(indexPointer) = m;
        indexPointer = indexPointer + 1;
    end
end

%% 3. Generate bits stream
bitsPerSymbol = log2(QAM_order);
tx_bits = randi([0,1], M_tot*bitsPerSymbol, 1);
x_vector = qammod(tx_bits, QAM_order, 'InputType', 'bit', 'UnitAveragePower', true);

%% 4. x_grid mapping by bits stream
% x_grid(k+1,l+1) represents the mathematical symbol x[k,l].
x_grid = complex(zeros(Nscale, M_scale(end)));
symbolPointer = 1;
for k = 0:Nscale-1
    Mk = M_scale(k+1);
    symbolIndices = symbolPointer:symbolPointer+Mk-1;
    x_grid(k+1,1:Mk) = x_vector(symbolIndices).';
    symbolPointer = symbolPointer + Mk;
end

%% 5. Discrete inverse Mellin-Fourier transform, Eq. (38)
% Direct summation implementation:
%
% X[n,m] = q^(-n/2)/Nscale *
%          sum_k 1/M(k) *
%          sum_l x[k,l] exp(j2pi(ml/M(k)-nk/Nscale))
X_grid = complex(zeros(Nscale, M_scale(end)));
for n = 0:Nscale-1
    Mn = M_scale(n+1);
    for m = 0:Mn-1
        temp = 0;
        for k = 0:Nscale-1
            Mk = M_scale(k+1);
            for l = 0:Mk-1
                temp = temp + x_grid(k+1,l+1)* exp(1j*2*pi*(m*l/Mk - n*k/Nscale))/Mk;
            end
        end
        X_grid(n+1,m+1) = q^(-n/2)*temp/Nscale;
    end
end

%% 6. Convert X_grid into X_vector
X_vector = complex(zeros(M_tot,1));
symbolPointer = 1;
for n = 0:Nscale-1
    Mn = M_scale(n+1);
    symbolIndices = symbolPointer:symbolPointer+Mn-1;
    X_vector(symbolIndices) = X_grid(n+1,1:Mn).';
    symbolPointer = symbolPointer + Mn;
end

%% 7. Construct T_iMF
% T_iMF is not used for the forward transform.
% It is retained for:
%   1. checking the direct Eq. (38) implementation;
%   2. applying the inverse transform at the receiver.
T_iMF = complex(zeros(M_tot,M_tot));

for outIdx = 1:M_tot
    n = pair_n(outIdx);
    m = pair_m(outIdx);

    for inIdx = 1:M_tot
        k = pair_n(inIdx);
        l = pair_m(inIdx);
        Mk = M_scale(k+1);

        T_iMF(outIdx,inIdx) = q^(-n/2)/(Nscale*Mk)*exp(1j*2*pi*(m*l/Mk - n*k/Nscale));
    end
end

X_vector_matrix = T_iMF*x_vector;

loopMatrixError = norm(X_vector-X_vector_matrix)/max(norm(X_vector_matrix),eps);

fprintf('Transform diagnostics\n');
fprintf('rank(T_iMF) = %d / %d\n', rank(T_iMF), M_tot);
fprintf('cond(T_iMF) = %.3e\n', cond(T_iMF));
fprintf('loop-versus-matrix error = %.3e\n\n', loopMatrixError);

%% 8. Build sampled ODSS synthesis matrix, Eq. (40)/(73)
Ns = round(T*Fs);
t = (0:Ns-1).'/Fs;

G_tx = complex(zeros(Ns,M_tot));
for col = 1:M_tot
    n = pair_n(col);
    m = pair_m(col);

    tLocal = q^n*(t - m/(q^n*W));

    G_tx(:,col) = q^(n/2)*odssTransmitPulse(tLocal, q, W, T, pulseType);
end

Gram_tx = Ts*(G_tx'*G_tx);
G_rx = G_tx/Gram_tx;
biorthErr = norm(Ts*(G_rx'*G_tx)-eye(M_tot),'fro')/sqrt(M_tot);

fprintf('Waveform diagnostics\n');
fprintf('samples = %d\n', Ns);
fprintf('rank(G_tx) = %d / %d\n', rank(G_tx), M_tot);
fprintf('cond(Ts*G_tx^H*G_tx) = %.3e\n', cond(Gram_tx));
fprintf('biorthogonality error = %.3e\n\n', biorthErr);

%% 9. ODSS waveform modulation
s_tx = G_tx*X_vector;

%% 10. Noiseless closed-loop check
Y0_vector = Ts*(G_rx'*s_tx);

x0_soft = T_iMF\Y0_vector;

fprintf('Noiseless closed-loop\n');
fprintf('||Y0-X||/||X|| = %.3e\n', norm(Y0_vector-X_vector)/norm(X_vector));
fprintf('||x0-x||/||x|| = %.3e\n\n', norm(x0_soft-x_vector)/norm(x_vector));

%% 11. Ideal no-delay, no-scale-spread channel plus AWGN
%
% h(tau,alpha) = h0*delta(tau)*delta(alpha-1)
%
% r(t) = h0*s(t) + w(t)

h0 = 1;%0.8 * exp(1j*pi/5);   % Example constant complex channel gain

s_channel = h0 * s_tx;

% Define SNR using received signal power after channel gain
signalPower = mean(abs(s_channel).^2);
noisePower = signalPower / 10^(SNRdB/10);

w_time = sqrt(noisePower/2) .* (randn(Ns,1) + 1j*randn(Ns,1));

r = s_channel + w_time;

measuredSNRdB = 10*log10(mean(abs(s_channel).^2) / mean(abs(w_time).^2));

fprintf('Ideal constant-gain channel plus AWGN\n');
fprintf('h0 = %.6f %+.6fj\n', real(h0), imag(h0));
fprintf('|h0| = %.6f\n', abs(h0));
fprintf('signal power = %.6e\n', signalPower);
fprintf('noise power  = %.6e\n', noisePower);
fprintf('target SNR   = %.2f dB\n', SNRdB);
fprintf('measured SNR = %.2f dB\n\n', measuredSNRdB);

%% 12. ODSS demodulator using matched-filter function, Eq. (45)
% Each column G_rx(:,idx) is one delay-scale receive template.
%
% Yhat[n,m] = integral conj(g_rx,nm(t)) * r(t) dt
Yhat_vector = complex(zeros(M_tot,1));
W_vector = complex(zeros(M_tot,1));

for idx = 1:M_tot

    n = pair_n(idx);
    m = pair_m(idx);

    alpha_nm = q^n;
    tau_nm = m/(alpha_nm*W);

    % The canonical-dual receive template for this (n,m)
    template_nm = G_rx(:,idx);

    % Matched-filter output for the received signal
    Yhat_vector(idx) = matchedFilterSample(r, t, template_nm);

    % Matched-filter output for the time-domain noise
    W_vector(idx) = matchedFilterSample(w_time, t, template_nm);

end


%% Verify matched-filter implementation
Yhat_matrix = Ts * (G_rx' * r);
W_matrix = Ts * (G_rx' * w_time);

matchedFilterImplementationError = ...
    norm(Yhat_vector - Yhat_matrix) / ...
    max(norm(Yhat_matrix), eps);

noiseFilterImplementationError = ...
    norm(W_vector - W_matrix) / ...
    max(norm(W_matrix), eps);

%% Verify ideal constant-gain channel relation
% Yhat = h0*X + W
linearityErr = ...
    norm(Yhat_vector - (h0*X_vector + W_vector)) / ...
    max(norm(Yhat_vector), eps);

fprintf('Delay-scale extraction\n');
fprintf('matched-filter versus matrix error = %.3e\n', ...
    matchedFilterImplementationError);

fprintf('noise-filter versus matrix error   = %.3e\n', ...
    noiseFilterImplementationError);

fprintf('||Yhat-(h0*X+W)||/||Yhat||         = %.3e\n\n', ...
    linearityErr);

%% 13. Construct Eq. (65) channel and apply one-tap MMSE
%
% Eq. (65):
%
% H_nm[n,m] =
% sum_p h_p *
% A_{g_rx,g_tx}(
%     q^n*(m/W*(alpha_p-1)-alpha_p*tau_p),
%     1/alpha_p)

%% 13.1 Channel path parameters
% Current constant-gain, zero-delay, zero-scale-spread channel

h_path     = h0;
tau_path   = 0;
alpha_path = 1;

h_path     = h_path(:);
tau_path   = tau_path(:);
alpha_path = alpha_path(:);

Path_num = numel(h_path);

if numel(tau_path) ~= Path_num || ...
        numel(alpha_path) ~= Path_num
    error(['h_path, tau_path, and alpha_path ', ...
           'must have the same length.']);
end

if any(alpha_path <= 0)
    error('All alpha_path values must be positive.');
end

%% 13.2 Prototype transmit and receive pulses

% Prototype transmit pulse g_tx(t)
g_tx_base = odssTransmitPulse(t, q, W, T, pulseType);

% Prototype receive pulse function
%
% This must be the g_rx(t) used by the paper receiver.
% For a matched-pulse receiver:
g_rx_fun = @(u) odssTransmitPulse( u, q, W, T, pulseType);

%% 13.3 Evaluate Eq. (65)
A_origin_raw = delay_scale_ambiguity( g_tx_base, t, g_rx_fun, 0, 1);

if abs(A_origin_raw) < 1e-12
    error('A_{g_rx,g_tx}(0,1) is too small for normalization.');
end

H_diagonal = complex(zeros(M_tot,1));

for idx = 1:M_tot

    n = pair_n(idx);
    m = pair_m(idx);

    H_nm = 0;

    for p = 1:Path_num

        tauAmbiguity = q^n * ...
            ( ...
              (m/W)*(alpha_path(p)-1) ...
              - alpha_path(p)*tau_path(p) ...
            );

        scaleAmbiguity = 1/alpha_path(p);

        A_grx_gtx_raw = delay_scale_ambiguity( ...
            g_tx_base, ...
            t, ...
            g_rx_fun, ...
            tauAmbiguity, ...
            scaleAmbiguity);

        % Enforce A_{g_rx,g_tx}(0,1) = 1
        A_grx_gtx = A_grx_gtx_raw / A_origin_raw;

        H_nm = H_nm + ...
            h_path(p)*A_grx_gtx;
    end

    H_diagonal(idx) = H_nm;
end

H = diag(H_diagonal);
D = H;

%% 13.4 Construct H and D

% Eq. (65) gives H_nm[n,m], i.e. diagonal coefficients
H = diag(H_diagonal);

% D is the diagonal part of H
D = diag(diag(H));

%% 13.5 MMSE loading

Es_delayScale = mean(abs(X_vector).^2);
sigmaW2_delayScale = mean(abs(W_vector).^2);

mmseLoading = sigmaW2_delayScale / Es_delayScale;

%% 13.6 Matrix-form one-tap MMSE

Zhat_vector = D' * ((D*D' + mmseLoading*eye(M_tot)) \ Yhat_vector);

%% 13.7 Equivalent element-wise form

d = diag(D);

Zhat_vector_scalar = ...
    conj(d) ./ ...
    (abs(d).^2 + mmseLoading) .* ...
    Yhat_vector;

mmseImplementationError = ...
    norm(Zhat_vector-Zhat_vector_scalar) / ...
    max(norm(Zhat_vector),eps);

%% 13.8 Diagnostics
A_origin_normalized = ...
    A_origin_raw / A_origin_raw;

expectedConstantChannel = ...
    h0 * ones(M_tot,1);

eq65ConstantChannelError = ...
    norm(H_diagonal-expectedConstantChannel) / ...
    max(norm(expectedConstantChannel),eps);

diagonalModelError = ...
    norm(Yhat_vector - ...
    (D*X_vector + W_vector)) / ...
    max(norm(Yhat_vector),eps);

fprintf('Eq. (65) channel diagnostics\n');

fprintf('raw A_grx,gtx(0,1)        = %.6e %+.6ej\n', ...
    real(A_origin_raw), imag(A_origin_raw));

fprintf('normalized A(0,1)         = %.6e %+.6ej\n', ...
    real(A_origin_normalized), ...
    imag(A_origin_normalized));

fprintf('minimum |H_nm[n,m]|       = %.6e\n', ...
    min(abs(H_diagonal)));

fprintf('maximum |H_nm[n,m]|       = %.6e\n', ...
    max(abs(H_diagonal)));

fprintf('Eq.65 constant-channel check = %.3e\n', ...
    eq65ConstantChannelError);

fprintf('Eq.65 diagonal-model error   = %.3e\n', ...
    diagonalModelError);

%% 14. Inverse ODSS transform and QAM demodulation, Eq. (70)
x_soft = T_iMF\Zhat_vector;

rx_bits = qamdemod(x_soft, QAM_order, 'OutputType', 'bit', 'UnitAveragePower', true);

numBitErrors = nnz(rx_bits ~= tx_bits);
BER = numBitErrors/numel(tx_bits);

x_hat = qammod(rx_bits, QAM_order, 'InputType', 'bit', 'UnitAveragePower', true);

numSymbolErrors = nnz(x_hat ~= x_vector);
SER = numSymbolErrors/M_tot;

fprintf('Recovered symbols\n');
fprintf('||Zhat-X||/||X|| = %.3e\n', norm(Zhat_vector-X_vector)/norm(X_vector));
fprintf('||xsoft-x||/||x|| = %.3e\n', norm(x_soft-x_vector)/norm(x_vector));
fprintf('symbol errors = %d / %d\n', numSymbolErrors, M_tot);
fprintf('SER = %.6e\n', SER);
fprintf('bit errors = %d / %d\n', numBitErrors, numel(tx_bits));
fprintf('BER = %.6e\n', BER);

%% 15. Theoretical BER of Gray-coded QPSK over AWGN

bitsPerSymbol = log2(QAM_order);

EbN0dB = SNRdB - 10*log10(bitsPerSymbol);
EbN0_linear = 10^(EbN0dB/10);

BER_theory = qfunc(sqrt(2*EbN0_linear));

fprintf('\nBER comparison\n');
fprintf('SNR = Es/N0      = %.2f dB\n', SNRdB);
fprintf('equivalent Eb/N0 = %.2f dB\n', EbN0dB);
fprintf('simulated BER    = %.6e\n', BER);
fprintf('QPSK theory BER  = %.6e\n', BER_theory);

%% 16. Constellation
figure;

plot(real(x_soft), imag(x_soft), 'o');
hold on;

plot(real(x_vector), imag(x_vector), 'kx', 'LineWidth', 1.5);

grid on;
axis equal;

xlabel('In-phase');
ylabel('Quadrature');

title(sprintf('ODSS over pure AWGN, SNR = %.1f dB', SNRdB));

legend('Before slicing', 'Transmitted 4-QAM', 'Location', 'best');

%% Local function: basic chirplet and optional PHYDYAS window
function g_tx = odssTransmitPulse(tLocal, q, W, T, pulseType)

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
            window = window + 2*(-1)^kk*A(kk).*cos(2*pi*kk*tv/(K*T));
        end

    otherwise
        error('Unknown pulseType: %s', pulseType);
end

g_tx(valid) = window.*g0;

end



function A = matchedFilterSample(r, t, template)
%MATCHEDFILTERSAMPLE
% Compute one sampled matched-filter output:
%
%   A = integral conj(template(t)) * r(t) dt
%
% For uniformly sampled signals:
%
%   A approximately equals Ts * template^H * r

    r = r(:);
    t = t(:);
    template = template(:);

    if numel(r) ~= numel(t)
        error('r and t must have the same number of samples.');
    end

    if numel(template) ~= numel(r)
        error('template and r must have the same number of samples.');
    end

    if numel(t) < 2
        error('The time vector must contain at least two samples.');
    end

    TsLocal = t(2)-t(1);

    % Matched-filter inner product
    % MATLAB '' performs conjugate transpose.
    A = TsLocal*(template'*r);

end


function A = delay_scale_ambiguity( ...
    r, t, g_rx_fun, tau, alpha)
%DELAY_SCALE_AMBIGUITY
%
% Computes:
%
% A_{g_rx,r}(tau,alpha)
%   = integral conj(g_rx(alpha*(t-tau))) ...
%              * sqrt(alpha) * r(t) dt
%
% g_rx_fun must accept arbitrary time arguments.

    r = r(:);
    t = t(:);

    if numel(r) ~= numel(t)
        error('r and t must have the same number of samples.');
    end

    if numel(t) < 2
        error('t must contain at least two samples.');
    end

    if alpha <= 0
        error('alpha must be positive.');
    end

    % Construct the delay-scale receive template directly.
    % No interp1 is required because g_rx_fun can be evaluated
    % directly at arbitrary time arguments.
    template = sqrt(alpha) .* ...
        g_rx_fun(alpha .* (t - tau));

    template = template(:);

    % Reuse the existing matched-filter function
    A = matchedFilterSample(r, t, template);
end