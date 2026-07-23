%% ODSS_AWGN_SmokeTest.m
% Initial ODSS transmitter/receiver smoke test over an AWGN channel.
%
% Main paper relations used:
%   Eq. (38)/(39): X = T_iMF * x
%   Eq. (40):      s(t) = sum_{n,m} X[n,m] q^(n/2)
%                         g_tx(q^n(t-m/(q^n W)))
%   Eq. (45):      Y[n,m] is obtained by ambiguity/matched-filter sampling
%   Eq. (50):      in an ideal no-delay/no-scale channel, Y[n,m] = X[n,m] + W[n,m]
%   Eq. (70):      x_hat = slicer(T_iMF^{-1} Z_hat)
%
% This is deliberately a DEBUG version, not yet the final low-complexity
% PHYDYAS receiver. A discrete bi-orthogonal receive bank is constructed so
% that G_rx' * G_tx / Fs = I. Therefore the noiseless waveform chain should
% close to machine precision before a delay-scale channel is introduced.
%
% Recommended debug order:
%   1) Run with Nscale = 4 and SNRdB = Inf.
%   2) Run with Nscale = 4 and SNRdB = 20 dB.
%   3) Change Nscale to 7.
%   4) Replace the numerical dual receive bank by a practical matched
%      PHYDYAS/chirplet receiver.

clear; clc; close all;
rng(20260720);

%% 1. ODSS parameters
q          = 2;          % dyadic delay-scale tiling
Nscale     = 4;          % start small; paper example later uses n = 0,...,6
B          = 1280;       % occupied baseband bandwidth, Hz
oversample = 8;
Fs         = oversample * B;  % 10240 Hz
SNRdB      = 20;         % set Inf for a noiseless first test

nSet   = 0:Nscale-1;
Mscale = floor(q.^nSet + 1e-12);  % M(n) = floor(q^n)
Mtot   = sum(Mscale);

% Paper bandwidth allocation: sum_n q^n W = B.
W  = B * (q - 1) / (q^Nscale - 1);
T0 = 1 / W;              % simple minimum-duration pulse used in this smoke test

% Stacking convention: n/k is the outer index, m/l is the inner index.
offset = [0, cumsum(Mscale)];

fprintf('================ ODSS AWGN smoke test ================\n');
fprintf('q = %.3f, Nscale = %d, Mtot = %d\n', q, Nscale, Mtot);
fprintf('B = %.2f Hz, W = %.6f Hz, T0 = %.6f s, Fs = %.1f Hz\n', ...
    B, W, T0, Fs);
fprintf('M(n) = %s\n', mat2str(Mscale));

%% 2. Generate unit-power QPSK symbols x[k,l]
bitsTx = randi([0, 1], 2*Mtot, 1);
x = ((1 - 2*bitsTx(1:2:end)) + 1j*(1 - 2*bitsTx(2:2:end))) / sqrt(2);

%% 3. ODSS preprocessing transform: Eq. (38)/(39)
TiMF = buildODSSiMFMatrix(q, Mscale);
TiMFInv = TiMF \ eye(Mtot);  % precompute the exact numerical inverse, as in Eq. (70)

X = TiMF * x;

fprintf('\nTransform diagnostics\n');
fprintf('rank(T_iMF) = %d / %d\n', rank(TiMF), Mtot);
fprintf('cond(T_iMF) = %.3e\n', cond(TiMF));
fprintf('inverse check ||T^{-1}T-I||_F/sqrt(M) = %.3e\n', ...
    norm(TiMFInv*TiMF - eye(Mtot), 'fro')/sqrt(Mtot));

if cond(TiMF) > 1e3
    warning(['T_iMF is ill-conditioned for this grid. This is one reason ', ...
             'to start the AWGN debug with fewer scales.']);
end

%% 4. Continuous-time sampled ODSS modulator: Eq. (40)
% The normalized chirplet in Eq. (71) sweeps from 1/sqrt(q) to sqrt(q).
% The frequency scale below converts that normalized sweep into a physical
% base-pulse bandwidth W. For q = 2, the base chirplet sweeps W -> 2W.
freqScale = W * sqrt(q) / (q - 1);
f1 = freqScale / sqrt(q);
f2 = freqScale * sqrt(q);
chirpRate = (f2 - f1) / T0;

Ntx = ceil(T0 * Fs);
t = (0:Ntx-1).' / Fs;

Gtx = buildODSSWaveformMatrix(t, Fs, q, W, T0, f1, chirpRate, Mscale);
sTx = Gtx * X;

fprintf('\nWaveform diagnostics\n');
fprintf('Tx duration = %.6f s, samples = %d\n', Ntx/Fs, Ntx);
fprintf('rank(G_tx) = %d / %d\n', rank(Gtx), Mtot);
fprintf('signal average power = %.6e\n', mean(abs(sTx).^2));

%% 5. Construct a discrete bi-orthogonal receive bank
% Let D = integral G_tx^H G_tx dt ~= G_tx^H G_tx / Fs.
% Choose G_rx = G_tx D^{-1}; then G_rx^H G_tx / Fs = I.
% This is an AWGN debugging device. It isolates transform/indexing errors
% before replacing G_rx by a practical shifted/compressed pulse bank.
GramTx = (Gtx' * Gtx) / Fs;
Grx = Gtx / GramTx;
biorthError = norm((Grx' * Gtx)/Fs - eye(Mtot), 'fro') / sqrt(Mtot);

fprintf('cond(G_tx^H G_tx / Fs) = %.3e\n', cond(GramTx));
fprintf('bi-orthogonality error = %.3e\n', biorthError);

%% 6. Noiseless receiver check
% Eq. (45): every entry is a sampled cross-ambiguity/matched-filter output.
Y0 = (Grx' * sTx) / Fs;
x0Linear = TiMFInv * Y0;
x0Slice = qpskSlice(x0Linear);
bits0 = qpskDemap(x0Slice);

fprintf('\nNoiseless closed-loop check\n');
fprintf('||Y0-X||/||X|| = %.3e\n', norm(Y0-X)/max(norm(X), eps));
fprintf('||x0-x||/||x|| = %.3e\n', norm(x0Linear-x)/norm(x));
fprintf('noiseless BER = %.6g (%d/%d)\n', ...
    mean(bits0 ~= bitsTx), nnz(bits0 ~= bitsTx), numel(bitsTx));

%% 7. AWGN channel: r(t) = s(t) + w(t)
signalPower = mean(abs(sTx).^2);
if isinf(SNRdB)
    noiseVar = 0;
else
    noiseVar = signalPower / (10^(SNRdB/10));
end

w = sqrt(noiseVar/2) * (randn(size(sTx)) + 1j*randn(size(sTx)));
r = sTx + w;

measuredSNRdB = 10*log10(sum(abs(sTx).^2) / max(sum(abs(w).^2), realmin));

%% 8A. Paper-structured AWGN receiver
% In the identity channel, Eq. (50) gives Y = X + W. Therefore D = I and
% no channel equalization is needed; invert T_iMF and slice to QPSK.
Yhat = (Grx' * r) / Fs;
Zhat = Yhat;
xHatLinear = TiMFInv * Zhat;
xHat = qpskSlice(xHatLinear);
bitsHat = qpskDemap(xHat);

berPaper = mean(bitsHat ~= bitsTx);
relErrPaper = norm(xHatLinear - x) / norm(x);

%% 8B. Optional full end-to-end LMMSE debug receiver
% This is not the final low-complexity ODSS receiver. It is useful as a
% reference: if this works but 8A fails, inspect the receive bank / inverse
% transform rather than the transmitted waveform.
Aend = Gtx * TiMF;
xHatLMMSELinear = (Aend' * Aend + noiseVar*eye(Mtot)) \ (Aend' * r);
xHatLMMSE = qpskSlice(xHatLMMSELinear);
bitsHatLMMSE = qpskDemap(xHatLMMSE);

berLMMSE = mean(bitsHatLMMSE ~= bitsTx);
relErrLMMSE = norm(xHatLMMSELinear - x) / norm(x);

fprintf('\nAWGN results\n');
fprintf('target SNR = %.2f dB, measured SNR = %.2f dB\n', SNRdB, measuredSNRdB);
fprintf('paper-structured receiver: BER = %.6g (%d/%d), rel. error = %.3e\n', ...
    berPaper, nnz(bitsHat ~= bitsTx), numel(bitsTx), relErrPaper);
fprintf('full end-to-end LMMSE:     BER = %.6g (%d/%d), rel. error = %.3e\n', ...
    berLMMSE, nnz(bitsHatLMMSE ~= bitsTx), numel(bitsTx), relErrLMMSE);

%% 9. Useful plots
figure('Name','ODSS transmit waveform');
plot(t, real(sTx));
grid on;
xlabel('Time (s)');
ylabel('Real\{s(t)\}');
title('ODSS transmit waveform');

Nfft = 2^nextpow2(max(4096, 8*Ntx));
S = fft(sTx, Nfft);
f = (0:Nfft-1).' * Fs/Nfft;
SdB = 20*log10(abs(S)/max(abs(S)) + eps);

figure('Name','ODSS transmit spectrum');
plot(f, SdB);
grid on;
xlim([0, min(Fs/2, 1.25*B)]);
ylim([-80, 5]);
xlabel('Frequency (Hz)');
ylabel('Normalized magnitude (dB)');
title('ODSS transmit spectrum');

figure('Name','Recovered QPSK symbols');
plot(real(xHatLinear), imag(xHatLinear), 'o');
hold on;
plot(real(x), imag(x), 'x', 'LineWidth', 1.5);
grid on; axis equal;
xlabel('In-phase');
ylabel('Quadrature');
legend('Before slicing', 'Transmitted', 'Location', 'best');
title(sprintf('Paper-structured receiver, SNR = %.1f dB', SNRdB));

%% ======================== Local functions ==============================
function TiMF = buildODSSiMFMatrix(q, Mscale)
%BUILDODSSIMFMATRIX Matrix implementation of paper Eq. (38).

    N = numel(Mscale);
    Mtot = sum(Mscale);
    offset = [0, cumsum(Mscale)];
    TiMF = complex(zeros(Mtot, Mtot));

    for nIdx = 1:N
        n = nIdx - 1;
        for m = 0:Mscale(nIdx)-1
            row = offset(nIdx) + m + 1;

            for kIdx = 1:N
                k = kIdx - 1;
                Mk = Mscale(kIdx);

                for l = 0:Mk-1
                    col = offset(kIdx) + l + 1;
                    phase = exp(1j*2*pi*(m*l/Mk - n*k/N));
                    TiMF(row, col) = q^(-n/2) * phase / (N*Mk);
                end
            end
        end
    end
end

function Gtx = buildODSSWaveformMatrix(t, Fs, q, W, T0, f1, chirpRate, Mscale)
%BUILDODSSWAVEFORMMATRIX Sampled waveform matrix for paper Eq. (40).
% Each column is q^(n/2) g_tx(q^n(t-m/(q^n W))).

    N = numel(Mscale);
    Mtot = sum(Mscale);
    offset = [0, cumsum(Mscale)];
    Gtx = complex(zeros(numel(t), Mtot));

    for nIdx = 1:N
        n = nIdx - 1;
        alpha = q^n;

        for m = 0:Mscale(nIdx)-1
            colIdx = offset(nIdx) + m + 1;
            tau = m / (alpha*W);
            tLocal = alpha * (t - tau);
            valid = (tLocal >= 0) & (tLocal < T0);

            g0 = complex(zeros(size(t)));
            g0(valid) = exp(1j*2*pi*( ...
                f1*tLocal(valid) + 0.5*chirpRate*tLocal(valid).^2)) / sqrt(T0);

            basis = sqrt(alpha) * g0;

            % Correct small finite-sampling energy differences. This keeps
            % every sampled basis waveform at unit continuous-time energy.
            energy = sum(abs(basis).^2) / Fs;
            if energy <= eps
                error('A waveform basis column has zero energy. Increase Fs or check parameters.');
            end
            Gtx(:, colIdx) = basis / sqrt(energy);
        end
    end
end

function xSlice = qpskSlice(z)
%QPSKSLICE Nearest-neighbour slicer for unit-power QPSK.
    bI = real(z) < 0;
    bQ = imag(z) < 0;
    xSlice = ((1 - 2*bI) + 1j*(1 - 2*bQ)) / sqrt(2);
end

function bits = qpskDemap(x)
%QPSKDEMAP Inverse of the QPSK mapping used in the main script.
    bits = zeros(2*numel(x), 1);
    bits(1:2:end) = real(x) < 0;
    bits(2:2:end) = imag(x) < 0;
end
