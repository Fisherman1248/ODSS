%% ODSS Modulation
% This self-contained example implements a small Orthogonal Delay Scale
% Space (ODSS) communication link in a tutorial style similar to the
% MathWorks OTFS example.
%
% Mathematical source:
% A. K. P. and C. R. Murthy, "Orthogonal Delay Scale Space Modulation:
% A New Technique for Wideband Time-Varying Channels," IEEE TSP, 2022.
%
% Implemented paper equations:
%   (1)  wideband delay-scale channel
%   (38) inverse Mellin-Fourier ODSS transform
%   (39) vector form X = TiMF*x
%   (40) ODSS waveform synthesis
%   (45)-(46) sampled delay-scale ambiguity / matched-filter receiver
%   (68)-(70) delay-scale equalization and inverse transform
%   (71)-(73) chirplet pulse and q-adic ODSS basis functions
%
% IMPORTANT:
% 1. This is a readable reference implementation, not a reproduction of
%    every BER curve in the paper.
% 2. The paper writes symbol recovery as TiMF^{-1} in Eq. (70). This
%    example therefore uses a MATLAB linear solve (TiMF \ X); it does not
%    assume that TiMF' is the inverse.
% 3. Exact robust bi-orthogonality is impossible with practical finite
%    pulses. Full LMMSE/least-squares recovery is used as the correctness
%    reference. The one-tap receiver is also shown as an approximation.
% 4. This tutorial maps the normalized basic-chirplet interval to
%    f1 = W/(q-1), f2 = qW/(q-1), so that f2-f1 = W and successive
%    q-adic compressions occupy successively wider frequency intervals.
%    This physical-frequency convention is documented explicitly.

clear;
clc;
close all;

%% Simulation Setup
rng(7, "twister");

cfg.q = 2;                  % dyadic scale ratio
cfg.Nscale = 7;             % n = 0,...,6, as in the paper example
cfg.B = 1280;               % system bandwidth in Hz
cfg.W = cfg.B*(cfg.q-1)/(cfg.q^cfg.Nscale-1); % Eq. (77)
cfg.T = 1.9;                % ODSS pulse/block duration in seconds
cfg.oversampling = 8;       % paper-style waveform oversampling factor
cfg.Fs = cfg.oversampling*cfg.B; % internal waveform sample rate in Hz
cfg.pulseType = "phydyas";     % keep "rect" while validating the basis
cfg.SNRdB = 25;             % waveform sample SNR for the channel example

% Single-path delay-scale channel. Change these values independently to
% study pure delay, pure scale, or their combination.
chan.h = exp(1j*pi/7);
chan.tau = 2/cfg.Fs;        % two-sample propagation delay
chan.alpha = 1.001;         % paper-style small time-scale factor

fprintf("ODSS tutorial configuration\n");
fprintf("  q = %.4g, Nscale = %d, B = %.1f Hz, W = %.6f Hz\n", ...
    cfg.q, cfg.Nscale, cfg.B, cfg.W);
fprintf("  T = %.3f s, Fs = %.1f Hz, oversampling = %d, pulse = %s\n\n", ...
    cfg.T, cfg.Fs, cfg.oversampling, cfg.pulseType);

%% ODSS Delay-Scale Grid
% The published paper uses M(n) = ceil(q^n). For q = 2 this also equals
% floor(q^n), but ceil is used here to follow the typeset equation.
%
% Vectorization order:
%   (n,m) = (0,0), (1,0),...,(1,M(1)-1),...,(N-1,M(N-1)-1)

grid = odssIndexMap(cfg.q, cfg.Nscale, cfg.W);

fprintf("Scale-layer sizes M(n): ");
fprintf("%d ", grid.M);
fprintf("\nTotal symbols Mtot = %d\n\n", grid.Mtot);

figure("Name", "ODSS delay-scale grid");
scatter(1e3*grid.tau, grid.alpha, 32, grid.n, "filled");
grid on;
xlabel("Delay grid point m/(q^n W) (ms)");
ylabel("Scale grid point q^n");
title(sprintf("ODSS irregular delay-scale grid, M_{tot} = %d", grid.Mtot));
colorbar;

%% Data Generation
% BPSK is used so the example requires no Communications Toolbox.

x = 2*randi([0 1], grid.Mtot, 1) - 1;
assert(all(abs(x) == 1), "BPSK generation failed.");

%% Inverse Mellin-Fourier Transform
% Equation (38):
%
% X[n,m] = q^(-n/2)/N * sum_k {
%              (1/M(k)) * sum_l x[k,l] *
%              exp(j*2*pi*(m*l/M(k) - n*k/N))
%          }.
%
% The literal nested-sum implementation is compared with a matrix form.

Xdirect = odssInverseMFDirect(x, grid, cfg.q, cfg.Nscale);
TiMF = buildTiMF(grid, cfg.q, cfg.Nscale);
X = TiMF*x;

transformAgreement = norm(Xdirect-X)/max(norm(Xdirect), eps);
fprintf("Direct-vs-matrix inverse-MF error = %.3e\n", transformAgreement);
fprintf("cond(TiMF) = %.3e\n", cond(TiMF));
assert(transformAgreement < 1e-11, ...
    "Direct and matrix implementations of Eq. (38) disagree.");

% Equation (70) explicitly uses TiMF^{-1}; use a solve, not TiMF'.
xTransformCheck = TiMF\X;
transformInverseError = norm(xTransformCheck-x)/norm(x);
fprintf("Transform inversion error = %.3e\n\n", transformInverseError);
assert(transformInverseError < 1e-10, ...
    "TiMF solve did not recover the source symbols.");

%% ODSS Transmit Basis
% Equations (71)-(73):
%
% g0(t) = exp(j*2*pi*(f1*t + 0.5*kappa*t^2))
%
% g_nm(t) = q^(n/2) * gtx(q^n*(t - m/(q^n*W))).
%
% Each column of Gtx is one sampled ODSS basis waveform g_nm(t).

tTx = (0:round(cfg.T*cfg.Fs)).'/cfg.Fs;
[Gtx, basisInfo] = odssTransmitBasis(tTx, grid, cfg);
dt = 1/cfg.Fs;

basisEnergies = real(diag(dt*(Gtx'*Gtx)));
fprintf("Sampled basis energy range = [%.6f, %.6f]\n", ...
    min(basisEnergies), max(basisEnergies));

Gram = dt*(Gtx'*Gtx);
basisNorm = sqrt(max(real(diag(Gram)), eps));
GramNormalized = Gram./(basisNorm*basisNorm.');
offDiagOnly = GramNormalized-diag(diag(GramNormalized));
offDiagGramRatio = norm(offDiagOnly, "fro") / ...
                   norm(diag(diag(GramNormalized)), "fro");

offMask = ~eye(grid.Mtot);
offValues = GramNormalized(offMask);
maxCrossCorr = max(abs(offValues));
rmsCrossCorr = sqrt(mean(abs(offValues).^2));
maxCrossCorrdB = 20*log10(maxCrossCorr + eps);
rmsCrossCorrdB = 20*log10(rmsCrossCorr + eps);

Cworst = abs(GramNormalized);
Cworst(1:grid.Mtot+1:end) = 0;
[worstCorr, worstLinearIndex] = max(Cworst(:));
[worstRxIndex, worstTxIndex] = ind2sub(size(Cworst), worstLinearIndex);

fprintf("\nWorst-correlated basis pair:\n");
fprintf("  Basis 1: index=%d, n=%d, m=%d\n", ...
    worstRxIndex, grid.n(worstRxIndex), grid.m(worstRxIndex));
fprintf("  Basis 2: index=%d, n=%d, m=%d\n", ...
    worstTxIndex, grid.n(worstTxIndex), grid.m(worstTxIndex));
fprintf("  Correlation magnitude = %.6f (%.2f dB)\n", ...
    worstCorr, 20*log10(worstCorr + eps));

fprintf("Normalized basis off-diagonal Frobenius ratio = %.3e\n", ...
    offDiagGramRatio);
fprintf("Maximum normalized cross-correlation = %.2f dB\n", ...
    maxCrossCorrdB);
fprintf("RMS normalized cross-correlation = %.2f dB\n\n", ...
    rmsCrossCorrdB);

figure("Name", "ODSS basis diagnostics");
subplot(1,2,1);
imagesc(20*log10(abs(GramNormalized)+1e-12));
axis image;
colorbar;
caxis([-60 0]);
xlabel("Basis index");
ylabel("Basis index");
title("Normalized basis Gram matrix (dB)");

subplot(1,2,2);
hold on;
representative = find(grid.m == 0);
nfftPlot = 2^nextpow2(4*numel(tTx));
f = (0:nfftPlot-1)*(cfg.Fs/nfftPlot);
for ii = 1:numel(representative)
    spectrum = abs(fft(Gtx(:,representative(ii)), nfftPlot));
    spectrum = spectrum/max(spectrum+eps);
    plot(f, 20*log10(spectrum+1e-12), ...
        "DisplayName", sprintf("n=%d", grid.n(representative(ii))));
end
grid on;
xlim([0 cfg.B]);
ylim([-70 3]);
xlabel("Frequency (Hz)");
ylabel("Normalized magnitude (dB)");
title("Representative ODSS basis spectra");
legend("Location", "best");

%% ODSS Modulation
% Equation (40), in matrix form:
%
%   s = Gtx*X = Gtx*TiMF*x.

s = Gtx*X;
sourceToWaveform = Gtx*TiMF;
assert(norm(s-sourceToWaveform*x)/max(norm(s),eps) < 1e-12, ...
    "Waveform matrix relation is inconsistent.");

figure("Name", "ODSS transmit waveform");
plot(tTx, real(s), "DisplayName", "Real");
hold on;
plot(tTx, imag(s), "DisplayName", "Imaginary");
grid on;
xlabel("Time (s)");
ylabel("Amplitude");
title("ODSS transmitted waveform");
legend;

%% Identity-Channel Verification
% Exact robust bi-orthogonality is not expected for finite practical
% pulses. Therefore the identity-channel correctness reference uses the
% complete waveform matrix, not a diagonal approximation.

xHatIdentity = fullLMMSE(sourceToWaveform, s, 0);
berIdentity = mean((real(xHatIdentity) >= 0) ~= (x >= 0));
identityRelativeError = norm(xHatIdentity-x)/norm(x);

fprintf("Identity-channel full-receiver relative error = %.3e\n", ...
    identityRelativeError);
fprintf("Identity-channel BER = %.3g\n\n", berIdentity);

assert(berIdentity == 0, ...
    "Identity-channel/no-noise BER must be zero.");
assert(identityRelativeError < 1e-9, ...
    "Identity-channel source recovery error is too large.");

%% Stop here while validating the transmit basis
fprintf("\nBasis-validation stage complete.\n");
fprintf("Proceeding to channel and receiver validation.\n");
% return;  % Uncomment this line to stop before the channel section.

%% Single-Path Delay-Scale Channel
% Equation (1), specialized to one discrete path:
%
%   r_s(t) = h*sqrt(alpha)*s(alpha*(t-tau)).
%
% The interpolation is a sampled approximation. Samples outside the
% transmit interval are set to zero.

tRxEnd = chan.tau + cfg.T/chan.alpha + 4/cfg.Fs;
tRx = (0:ceil(tRxEnd*cfg.Fs)).'/cfg.Fs;

% Apply the same physical channel to every ODSS basis column. Gch maps
% delay-scale coefficients X to the received waveform.
Gch = applyDelayScaleToBasis(Gtx, tTx, tRx, chan);
rs = Gch*X;

fprintf("Channel: |h| = %.3f, angle(h) = %.3f rad, tau = %.3f ms, alpha = %.6f\n", ...
    abs(chan.h), angle(chan.h), 1e3*chan.tau, chan.alpha);

%% AWGN
% SNR here is defined as:
%
%   mean(|r_s|^2) / E{|w|^2}
%
% at the sampled receiver waveform. It is not labeled Eb/N0.

signalPower = mean(abs(rs).^2);
noiseVar = signalPower/10^(cfg.SNRdB/10);
w = sqrt(noiseVar/2)*(randn(size(rs)) + 1j*randn(size(rs)));
r = rs + w;

measuredSNRdB = 10*log10(mean(abs(rs).^2)/mean(abs(w).^2));
fprintf("Requested waveform SNR = %.2f dB, measured = %.2f dB\n\n", ...
    cfg.SNRdB, measuredSNRdB);

%% Full LMMSE Equalization
% The complete source-to-received-waveform matrix is:
%
%   A = Gch*TiMF.
%
% The full waveform-domain LMMSE estimate is:
%
%   xHat = (A^H*A + sigma^2*I)^(-1) A^H*r.

Afull = Gch*TiMF;
xHatFull = fullLMMSE(Afull, r, noiseVar);
berFull = mean((real(xHatFull) >= 0) ~= (x >= 0));

%% ODSS Receiver and One-Tap Equalization
% Equations (45)-(46) sample the scale-delay cross-ambiguity function.
% With grx = gtx, its sampled approximation is a bank of inner products:
%
%   Yhat = dt*Grx^H*r.
%
% Dfull below is the sampled delay-scale effective channel. Equation (69)
% is applied to diag(Dfull) to demonstrate the paper's low-complexity
% one-tap approximation.

[Grx, ~] = odssTransmitBasis(tRx, grid, cfg);
Yhat = dt*(Grx'*r);
Dfull = dt*(Grx'*Gch);

%% Full Matched-Filter-Domain Equalization
% Matched-filter input-output relation:
%
%   Yhat = Dfull*TiMF*x + w_MF.
%
% This receiver retains the complete Dfull matrix and therefore serves as
% the delay-scale-domain correctness reference. The one-tap receiver below
% keeps only diag(Dfull).

Amatched = Dfull*TiMF;
receiveBasisEnergy = real(diag(dt*(Grx'*Grx)));
noiseVarMatchedEach = noiseVar*dt*receiveBasisEnergy;
noiseVarMatched = mean(noiseVarMatchedEach);

if cfg.SNRdB >= 80
    xHatMatchedFull = Amatched\Yhat;
else
    xHatMatchedFull = ...
        (Amatched'*Amatched + noiseVarMatched*eye(grid.Mtot)) \ ...
        (Amatched'*Yhat);
end

berMatchedFull = mean((real(xHatMatchedFull) >= 0) ~= (x >= 0));
matchedFullError = norm(xHatMatchedFull-x)/norm(x);

fprintf("Full matched-filter-domain relative error = %.3e\n", ...
    matchedFullError);
fprintf("Full matched-filter-domain BER = %.3e\n", ...
    berMatchedFull);

[XHatOneTap, noiseVarDelayScale] = oneTapLMMSE( ...
    Yhat, Dfull, noiseVar, Grx, dt);

% Equation (70): xHat = slice(TiMF^{-1}*Zhat).
xHatOneTapSoft = TiMF\XHatOneTap;
berOneTap = mean((real(xHatOneTapSoft) >= 0) ~= (x >= 0));

Ddiag = diag(diag(Dfull));
offDiagonalEnergyRatio = norm(Dfull-Ddiag, "fro")^2 / ...
                         max(norm(Ddiag, "fro")^2, eps);

fprintf("Delay-scale effective-channel off-diagonal energy ratio = %.3e\n", ...
    offDiagonalEnergyRatio);
fprintf("Approximate matched-filter noise variance range = [%.3e, %.3e]\n", ...
    min(noiseVarDelayScale), max(noiseVarDelayScale));
fprintf("Full LMMSE BER = %.3e\n", berFull);
fprintf("One-tap delay-scale BER = %.3e\n\n", berOneTap);

if berOneTap > 0.05
    warning(['The diagonal one-tap approximation has a high BER in this ' ...
        'configuration. The full waveform and full matched-filter-domain ' ...
        'receivers remain the correctness references.']);
end

figure("Name", "ODSS received-channel diagnostics");
subplot(1,2,1);
Dscale = max(abs(Dfull(:)));
imagesc(20*log10(abs(Dfull)/max(Dscale,eps)+1e-12));
axis image;
colorbar;
caxis([-60 0]);
xlabel("Transmit delay-scale index");
ylabel("Receive delay-scale index");
title("Effective delay-scale channel D (dB)");

subplot(1,2,2);
plot(real(x), imag(x), "ko", "DisplayName", "Transmitted BPSK");
hold on;
plot(real(xHatFull), imag(xHatFull), "b+", ...
    "DisplayName", "Waveform full LMMSE");
plot(real(xHatMatchedFull), imag(xHatMatchedFull), "gd", ...
    "DisplayName", "MF-domain full");
plot(real(xHatOneTapSoft), imag(xHatOneTapSoft), "rx", ...
    "DisplayName", "One-tap");
grid on;
axis equal;
xlabel("In-phase");
ylabel("Quadrature");
title("Detected source symbols");
legend("Location", "best");

%% Summary
fprintf("------------------------------------------------------------\n");
fprintf("ODSS example summary\n");
fprintf("  Mtot                                  : %d\n", grid.Mtot);
fprintf("  Eq. (38) direct/matrix error          : %.3e\n", transformAgreement);
fprintf("  TiMF inversion error                  : %.3e\n", transformInverseError);
fprintf("  Identity-channel BER                  : %.3e\n", berIdentity);
fprintf("  Channel waveform full-LMMSE BER       : %.3e\n", berFull);
fprintf("  Channel MF-domain full BER            : %.3e\n", berMatchedFull);
fprintf("  Channel MF-domain relative error      : %.3e\n", matchedFullError);
fprintf("  Channel one-tap BER                   : %.3e\n", berOneTap);
fprintf("  Basis maximum cross-correlation       : %.2f dB\n", ...
    maxCrossCorrdB);
fprintf("  Basis RMS cross-correlation           : %.2f dB\n", ...
    rmsCrossCorrdB);
fprintf("  Delay-scale off-diagonal energy ratio : %.3e\n", ...
    offDiagonalEnergyRatio);
fprintf("------------------------------------------------------------\n");

%% Local Support Functions

function grid = odssIndexMap(q, Nscale, W)
%ODSSINDEXMAP Construct the irregular ODSS scale-delay grid.
%
% Mathematical indices are zero-based. MATLAB vector indices are one-based.
% M(n) follows the published typeset definition M(n) = ceil(q^n).

    arguments
        q (1,1) double {mustBeGreaterThanOrEqual(q,1)}
        Nscale (1,1) double {mustBeInteger,mustBePositive}
        W (1,1) double {mustBePositive}
    end

    M = ceil(q.^(0:Nscale-1));
    Mtot = sum(M);

    nIndex = zeros(Mtot,1);
    mIndex = zeros(Mtot,1);
    offset = [0 cumsum(M(1:end-1))];

    cursor = 1;
    for n = 0:Nscale-1
        for m = 0:M(n+1)-1
            nIndex(cursor) = n;
            mIndex(cursor) = m;
            cursor = cursor + 1;
        end
    end

    grid.q = q;
    grid.Nscale = Nscale;
    grid.W = W;
    grid.M = M;
    grid.Mtot = Mtot;
    grid.offset = offset;
    grid.n = nIndex;
    grid.m = mIndex;
    grid.alpha = q.^nIndex;
    grid.tau = mIndex./(q.^nIndex*W); % paper Eqs. (40), (42), and (73)
end

function X = odssInverseMFDirect(x, grid, q, Nscale)
%ODSSINVERSEMFDIRECT Literal nested-sum implementation of paper Eq. (38).
%
% x and X are Mtot-by-1 vectors using the grid's scale-major ordering.

    validateattributes(x, {'double'}, {'column','numel',grid.Mtot});
    X = complex(zeros(grid.Mtot,1));

    for outIndex = 1:grid.Mtot
        n = grid.n(outIndex);
        m = grid.m(outIndex);
        accumulator = 0;

        for inIndex = 1:grid.Mtot
            k = grid.n(inIndex);
            l = grid.m(inIndex);
            Mk = grid.M(k+1);

            phase = exp(1j*2*pi*(m*l/Mk - n*k/Nscale));
            accumulator = accumulator + x(inIndex)*phase/Mk;
        end

        X(outIndex) = q^(-n/2)*accumulator/Nscale;
    end
end

function TiMF = buildTiMF(grid, q, Nscale)
%BUILDTIMF Construct the Mtot-by-Mtot matrix from paper Eq. (38).

    TiMF = complex(zeros(grid.Mtot, grid.Mtot));

    for outIndex = 1:grid.Mtot
        n = grid.n(outIndex);
        m = grid.m(outIndex);

        for inIndex = 1:grid.Mtot
            k = grid.n(inIndex);
            l = grid.m(inIndex);
            Mk = grid.M(k+1);

            TiMF(outIndex,inIndex) = ...
                q^(-n/2)/Nscale/Mk * ...
                exp(1j*2*pi*(m*l/Mk - n*k/Nscale));
        end
    end
end

function [G, info] = odssTransmitBasis(t, grid, cfg)
%ODSSTRANSMITBASIS Sample the ODSS basis in paper Eqs. (71)-(73).
%
% G is length(t)-by-Mtot. Column i corresponds to
% (n,m) = (grid.n(i), grid.m(i)).
%
% Physical-frequency convention used here:
%   f1 = W/(q-1), f2 = q*W/(q-1), hence f2-f1 = W.
% Successive q-adic compressions then expand the occupied bandwidth.

    t = t(:);
    q = cfg.q;
    W = cfg.W;
    T = cfg.T;
    Fs = cfg.Fs;

    f1 = W/(q-1);
    f2 = q*W/(q-1);
    kappa = (f2-f1)/T;

    % Normalize the unscaled base pulse to unit sampled energy.
    tBase = (0:round(T*Fs)).'/Fs;
    basePulse = localBasePulse(tBase, T, f1, kappa, cfg.pulseType);
    baseEnergy = sum(abs(basePulse).^2)/Fs;
    baseNorm = sqrt(baseEnergy);

    G = complex(zeros(numel(t), grid.Mtot));

    for index = 1:grid.Mtot
        n = grid.n(index);
        m = grid.m(index);

        delay = m/(q^n*W);
        u = q^n*(t-delay);
        support = (u >= 0) & (u <= T);

        column = complex(zeros(size(t)));
        column(support) = q^(n/2) * ...
            localBasePulse(u(support), T, f1, kappa, cfg.pulseType) / ...
            baseNorm;

        G(:,index) = column;
    end

    info.f1 = f1;
    info.f2 = f2;
    info.kappa = kappa;
    info.baseEnergyBeforeNormalization = baseEnergy;
end

function pulse = localBasePulse(t, T, f1, kappa, pulseType)
%LOCALBASEPULSE Construct rectangular or PHYDYAS-windowed chirplet.

    chirplet = exp(1j*2*pi*(f1*t + 0.5*kappa*t.^2));

    switch lower(char(pulseType))
        case "rect"
            window = ones(size(t));

        case "phydyas"
            % Paper Eq. (72), overlap factor K = 3.
            K = 3;
            A1 = 0.91143783;
            A2 = 0.41143783;
            window = 1 ...
                + 2*(-1)^1*A1*cos(2*pi*1*t/(K*T)) ...
                + 2*(-1)^2*A2*cos(2*pi*2*t/(K*T));

        otherwise
            error("Unsupported pulseType '%s'. Use 'rect' or 'phydyas'.", ...
                pulseType);
    end

    pulse = window.*chirplet;
end

function Gout = applyDelayScaleToBasis(Gin, tIn, tOut, chan)
%APPLYDELAYSCALETOBASIS Apply a discrete-path wideband channel to columns.
%
% Implements the one-path form of paper Eq. (1):
%   y(t) = h*sqrt(alpha)*x(alpha*(t-tau)).
%
% Linear interpolation is used; values outside tIn are zero.

    tIn = tIn(:);
    tOut = tOut(:);

    validateattributes(chan.alpha, {'double'}, {'scalar','positive'});
    validateattributes(chan.tau, {'double'}, {'scalar','nonnegative'});

    queryTime = chan.alpha*(tOut-chan.tau);
    Gout = complex(zeros(numel(tOut), size(Gin,2)));

    for column = 1:size(Gin,2)
        scaled = interp1(tIn, Gin(:,column), queryTime, "linear", 0);
        Gout(:,column) = chan.h*sqrt(chan.alpha)*scaled;
    end
end

function xHat = fullLMMSE(A, y, noiseVar)
%FULLLMMSE Full waveform-domain linear MMSE estimate.
%
% y = A*x + w, E[w*w^H] = noiseVar*I.

    if noiseVar < 0
        error("noiseVar must be nonnegative.");
    end

    if noiseVar == 0
        % More stable than forming normal equations in the noiseless case.
        xHat = A\y;
    else
        numSymbols = size(A,2);
        xHat = (A'*A + noiseVar*eye(numSymbols))\(A'*y);
    end
end

function [XHat, noiseVarDS] = oneTapLMMSE(Yhat, Dfull, noiseVar, Grx, dt)
%ONETAPLMMSE Diagonal delay-scale approximation of paper Eq. (69).
%
% Dfull is not forced to be diagonal. Only diag(Dfull) is retained.
% The post-matched-filter noise variance is approximated separately for
% every basis waveform.

    d = diag(Dfull);

    receiveBasisEnergy = real(diag(dt*(Grx'*Grx)));
    noiseVarDS = noiseVar*dt*receiveBasisEnergy;

    denominator = abs(d).^2 + noiseVarDS;
    denominator = max(denominator, eps);

    XHat = conj(d).*Yhat./denominator;
end
