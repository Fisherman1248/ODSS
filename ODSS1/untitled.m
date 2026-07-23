% ODSS code, generated from offical OTFS matlab example.
clc; clear; close all;

%% Simulation Set up
q = 2;
Nscale = 7;             % n = 0,...,6, as in the paper example
B = 1280;               % system bandwidth in Hz, fig. 8 illustrate
W = B*(q-1)/(q^Nscale-1); % Transmit filter bandwidth, Eq. (77)
T = 1.9;                % ODSS pulse/block duration in seconds
TimeScaleN = 7;
oversampling = 8; % fig. 9 illustrates
pulseType = 'phydyas_2';
QAM_order = 4;

SNRdB = 100;

% Waveform sampling rate
Fs = oversampling * B;
Ts = 1/Fs;

% Channel Simulation set up, follow fig. 9
Path_num = 1;
tau_max = 0;
alpha_max = 1.00;

tau_p = tau_max * rand(Path_num, 1); % Uniform delay: U(0, tau_max)

alpha_p = 1/alpha_max + (alpha_max - 1/alpha_max) * rand(Path_num, 1); % Uniform scale: U(1/alpha_max, alpha_max)

%% x_grid set up
n_set = 0 : Nscale-1;
M_scale = floor(q.^n_set); % The number of points in each scale layer
M_tot = sum(M_scale); % The total number of scale points

%% Generate bits stream
bitsPerSymbol = log2(QAM_order);    % 2 bits per 4-QAM symbol
% One ODSS block contains M_tot QAM symbols
tx_bits = randi([0, 1], M_tot * bitsPerSymbol, 1);
% Every two consecutive bits are mapped to one 4-QAM symbol
x_vector = qammod(tx_bits, QAM_order, 'InputType', 'bit', 'UnitAveragePower', true);

%% x_grid mapping by bits stream
% x(k+1,l+1) represents the mathematical symbol x[k,l]
x_grid = complex(zeros(Nscale, M_scale(end)));

symbolPointer = 1;

for k = 0:Nscale-1
    Mk = M_scale(k+1);

    symbolIndices = symbolPointer : symbolPointer + Mk - 1;

    % Valid grid points: l = 0,...,M(k)-1
    x_grid(k+1, 1:Mk) = x_vector(symbolIndices).';

    symbolPointer = symbolPointer + Mk;
end

%% X_gird (eq. 38)       X[n,m] <-- x[k,l]
X_grid = complex(zeros(Nscale, M_scale(end)));

% traverse for n in X[n,m]
for n = 0:Nscale-1
    Mn = M_scale(n+1);
    % traverse for m in X[n,m]
    for m = 0:Mn-1
        temp = 0;   
        % eq. 38 outer summation
        for k = 0:Nscale-1
            Mk = M_scale(k+1);
            l = 0:Mk-1;
            phase = exp(1j*2*pi * (m*l/Mk - n*k/Nscale));
            temp = temp + sum(x_grid(k+1, 1:Mk) .* phase) / Mk; % inner summation
        end
        X_grid(n+1, m+1) = q^(-n/2) * temp / Nscale; % scaling
    end
end

%% Convert X[n,m] to vector before waveform generation
X_vector = complex(zeros(M_tot,1));

symbolPointer = 1;

for n = 0:Nscale-1
    Mn = M_scale(n+1);

    X_vector(symbolPointer:symbolPointer+Mn-1) = ...
        X_grid(n+1,1:Mn).';

    symbolPointer = symbolPointer + Mn;
end

%% Generate sampled ODSS transmit basis and waveform
t = (0:Ts:T-Ts).';       % column vector

G_tx = complex(zeros(length(t), M_tot));

symbolPointer = 1;

for n = 0:Nscale-1
    Mn = M_scale(n+1);

    for m = 0:Mn-1
        tLocal = q^n * (t - m/(q^n*W));

        g_tx_nm = odssTransmitPulse( ...
            tLocal, q, W, T, pulseType);

        s_mn = q^(n/2) * g_tx_nm;

        G_tx(:,symbolPointer) = s_mn(:);

        symbolPointer = symbolPointer + 1;
    end
end

s_tx = G_tx * X_vector;

fprintf('Waveform sample rate: %.1f Hz\n', Fs);
fprintf('ODSS block duration: %.6f s\n', T);
fprintf('Number of samples: %d\n', length(s_tx));

%% Check individual ODSS scale spectra, similar to Fig. 8
Nfft = 2^nextpow2(length(t)) * 4;
f = (0:Nfft/2-1) * Fs/Nfft;
figure;
hold on;
for n = 0:Nscale-1
    % Use m = 0; time shifts only change spectral phase,
    % not the magnitude spectrum
    tLocal = q^n * t;
    g_tx_n = odssTransmitPulse(tLocal, q, W, T, pulseType);
    s_n = q^(n/2) * g_tx_n;
    S_n = fft(s_n, Nfft);
    S_n = S_n(1:Nfft/2);
    S_n_dB = 20*log10(abs(S_n)/max(abs(S_n)) + eps);
    plot(f, S_n_dB, 'DisplayName', sprintf('n = %d', n));
end
xlabel('Frequency (Hz)');
ylabel('Normalized magnitude (dB)');
title('Individual ODSS Scale Spectra');
xlim([0 B]);
ylim([-50 0]);
grid on;
legend('Location', 'best');


%% Pure AWGN channel
rs = s_tx;

Ps = mean(abs(rs).^2);
Pn = Ps * 10^(-SNRdB/10);

w = sqrt(Pn/2) * ...
    (randn(size(rs)) + 1j*randn(size(rs)));

r = rs + w;

measuredSNR = 10*log10( ...
    mean(abs(rs).^2) / mean(abs(w).^2));

fprintf('\nPure AWGN channel:\n');
fprintf('Signal power: %.6e\n', Ps);
fprintf('Noise power: %.6e\n', Pn);
fprintf('Target SNR: %.2f dB\n', SNRdB);
fprintf('Measured SNR: %.2f dB\n', measuredSNR);

%% Display channel information
fprintf('\nDelay-scale channel:\n');
fprintf('Number of paths: %d\n', Path_num);
fprintf('Receiver duration: %.6f s\n', tRx(end));
fprintf('Number of received samples: %d\n', length(r));
fprintf('Signal power: %.6e\n', Ps);
fprintf('Noise power: %.6e\n', Pn);

if ~isinf(SNRdB)
    measuredSNR = 10*log10( ...
        mean(abs(rs).^2) / mean(abs(w).^2));

    fprintf('Target SNR: %.2f dB\n', SNRdB);
    fprintf('Measured SNR: %.2f dB\n', measuredSNR);
end


%% Construct sampled biorthogonal receive basis
Gram = Ts * (G_tx' * G_tx);

fprintf('\nWaveform diagnostics:\n');
fprintf('rank(G_tx) = %d / %d\n', rank(G_tx), M_tot);
fprintf('cond(Gram) = %.3e\n', cond(Gram));

if rcond(Gram) < 1e-13
    error('G_tx Gram matrix is numerically singular.');
end

G_rx = G_tx / Gram;

biorthError = norm( ...
    Ts*(G_rx' * G_tx) - eye(M_tot), ...
    'fro') / sqrt(M_tot);

fprintf('Biorthogonality error = %.3e\n', ...
    biorthError);


%% Noiseless waveform closed-loop test
Y0_vector = Ts * G_rx' * s_tx;

Y0Error = norm(Y0_vector-X_vector) ...
    / max(norm(X_vector),eps);

fprintf('Noiseless ||Y0-X||/||X|| = %.3e\n', ...
    Y0Error);

%% AWGN ODSS receiver
Y_vector = Ts * G_rx' * r;

% Pure AWGN: equivalent channel is identity
Zhat_vector = Y_vector;


%% 9. QAM constellation slicing, Eq. (70)
rx_bits = qamdemod(x_hat_vector, QAM_order, 'OutputType', 'bit', 'UnitAveragePower', true);
rx_bits = rx_bits(:);

numBitErrors = sum(rx_bits ~= tx_bits);
BER = numBitErrors / numel(tx_bits);

fprintf('\nReceiver results:\n');
fprintf('Bit errors: %d / %d\n', numBitErrors, numel(tx_bits));
fprintf('BER: %.6e\n', BER);


%% 10. Plot transmitted and recovered QAM symbols
figure;
subplot(1,2,1);
scatter(real(x_vector), imag(x_vector), 25, 'filled');
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('Transmitted QAM Symbols');

subplot(1,2,2);
scatter(real(x_hat_vector), imag(x_hat_vector), 25, 'filled');
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('Recovered QAM Symbols');