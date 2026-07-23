%% odss_awgn_minimal.m
% Minimal ODSS baseband transceiver over AWGN
% No Communications Toolbox required.
%
% Chain:
% bits -> QPSK -> inverse Mellin-Fourier ODSS transform
% -> chirplet synthesis -> AWGN
% -> biorthogonal matched-filter bank
% -> inverse ODSS transform -> QPSK decision

clear; clc; close all;
rng(1);

%% 1. ODSS parameters
q       = 2;
Nscale  = 4;
nSet    = 0:Nscale-1;
Mscale  = floor(q.^nSet);
Mtot    = sum(Mscale);

B       = 1280;                         % nominal bandwidth [Hz]
W       = B*(q-1)/(q^Nscale-1);         % delay-grid parameter [Hz]
T       = 1/W;                          % basic pulse duration [s]
Fs      = 8*B;                          % sampling frequency [Hz]
Ns      = round(T*Fs);
t       = (0:Ns-1).'/Fs;

SNRdB      = 20;
numFrames  = 2000;

fprintf('q = %.1f, Nscale = %d, Mtot = %d\n', q, Nscale, Mtot);
fprintf('W = %.6f Hz, T = %.6e s, Fs = %.1f Hz, Ns = %d\n', ...
    W, T, Fs, Ns);

%% 2. Build input and output index maps
% Input vector x stacks x[k,l].
% Output vector X stacks X[n,m].

kMap = zeros(Mtot,1);
lMap = zeros(Mtot,1);
p = 1;

for k = 0:Nscale-1
    for l = 0:Mscale(k+1)-1
        kMap(p) = k;
        lMap(p) = l;
        p = p + 1;
    end
end

nMap = zeros(Mtot,1);
mMap = zeros(Mtot,1);
p = 1;

for n = 0:Nscale-1
    for m = 0:Mscale(n+1)-1
        nMap(p) = n;
        mMap(p) = m;
        p = p + 1;
    end
end

%% 3. Construct the inverse Mellin-Fourier ODSS transform
% Paper Eq. (121):
%
% X[n,m] = q^(-n/2)/N *
%          sum_k { 1/M(k) * sum_l x[k,l]
%          * exp(j*2*pi*(m*l/M(k) - n*k/N)) }

T_iMF = complex(zeros(Mtot,Mtot));

for outIdx = 1:Mtot
    n = nMap(outIdx);
    m = mMap(outIdx);

    for inIdx = 1:Mtot
        k = kMap(inIdx);
        l = lMap(inIdx);

        T_iMF(outIdx,inIdx) = ...
            q^(-n/2) / (Nscale*Mscale(k+1)) * ...
            exp(1j*2*pi*(m*l/Mscale(k+1) - n*k/Nscale));
    end
end

fprintf('rank(T_iMF) = %d / %d\n', rank(T_iMF), Mtot);
fprintf('cond(T_iMF) = %.3e\n', cond(T_iMF));

%% 4. Construct sampled ODSS chirplet synthesis matrix
% Basic rectangular-window chirplet used only for an AWGN smoke test.
% A discrete dual receive basis is constructed below so that
% G_rx^H * G_tx / Fs = I.

f1     = W/sqrt(q);
f2     = W*sqrt(q);
kappa  = (f2-f1)/T;

G_tx = complex(zeros(Ns,Mtot));

for col = 1:Mtot
    n = nMap(col);
    m = mMap(col);

    % Argument of g_tx(q^n(t-m/(q^n W)))
    u = q^n .* (t - m/(q^n*W));
    valid = (u >= 0) & (u < T);

    pulse = complex(zeros(Ns,1));

    % sqrt(W) gives the basic pulse approximately unit energy.
    pulse(valid) = sqrt(W)*q^(n/2) .* ...
        exp(1j*2*pi*(f1*u(valid) + 0.5*kappa*u(valid).^2));

    G_tx(:,col) = pulse;
end

Gram = (G_tx' * G_tx)/Fs;

if rcond(Gram) < 1e-12
    error('Waveform Gram matrix is numerically singular.');
end

% Canonical sampled biorthogonal receive bank.
G_rx = G_tx / Gram;

biorthError = norm((G_rx' * G_tx)/Fs - eye(Mtot),'fro')/sqrt(Mtot);

fprintf('rank(G_tx) = %d / %d\n', rank(G_tx), Mtot);
fprintf('cond(G_tx^H G_tx/Fs) = %.3e\n', cond(Gram));
fprintf('biorthogonality error = %.3e\n', biorthError);

%% 5. Noiseless closed-loop test
bits0 = randi([0 1],2*Mtot,1);

x0 = ((1-2*bits0(1:2:end)) + ...
    1j*(1-2*bits0(2:2:end)))/sqrt(2);

X0   = T_iMF*x0;
s0   = G_tx*X0;
Y0   = (G_rx' * s0)/Fs;
x0Hat = T_iMF \ Y0;

noiselessError = norm(x0Hat-x0)/norm(x0);

fprintf('noiseless relative error = %.3e\n', noiselessError);

%% 6. AWGN BER simulation
bitErrors = 0;
totalBits = 0;

for frame = 1:numFrames

    % QPSK mapping:
    % bit 0 -> +1 and bit 1 -> -1 on each real dimension.
    bits = randi([0 1],2*Mtot,1);

    x = ((1-2*bits(1:2:end)) + ...
        1j*(1-2*bits(2:2:end)))/sqrt(2);

    % ODSS preprocessing and waveform synthesis
    X = T_iMF*x;
    s = G_tx*X;

    % Complex AWGN at the requested sampled-waveform SNR
    signalPower = mean(abs(s).^2);
    noisePower  = signalPower / 10^(SNRdB/10);

    w = sqrt(noisePower/2) * ...
        (randn(Ns,1) + 1j*randn(Ns,1));

    r = s + w;

    % ODSS receiver
    YHat = (G_rx' * r)/Fs;
    xHatSoft = T_iMF \ YHat;

    % Hard QPSK decision
    bitsHat = zeros(2*Mtot,1);
    bitsHat(1:2:end) = real(xHatSoft) < 0;
    bitsHat(2:2:end) = imag(xHatSoft) < 0;

    bitErrors = bitErrors + sum(bitsHat ~= bits);
    totalBits = totalBits + numel(bits);
end

BER = bitErrors/totalBits;

fprintf('\nSNR = %.1f dB\n', SNRdB);
fprintf('frames = %d\n', numFrames);
fprintf('bit errors = %d / %d\n', bitErrors, totalBits);
fprintf('BER = %.6e\n', BER);

%% 7. Optional constellation of the last received block
figure;
plot(real(xHatSoft),imag(xHatSoft),'o','MarkerSize',6);
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title(sprintf('Recovered QPSK symbols, SNR = %.1f dB',SNRdB));
