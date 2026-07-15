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
pulseType = 'phydyas';
QAM_order = 4;

% Waveform sampling rate
Fs = oversampling * B;
Ts = 1/Fs;

% Channel Simulation set up, follow fig. 9
Path_num = 20;
tau_max = 0.01;
alpha_max = 1.001;

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


%% Generate ODSS transmitted waveform, Eq. (73)

% Time axis for one ODSS block
t = 0:Ts:(T-Ts);

% Initialize transmitted waveform
s_tx = complex(zeros(size(t)));

for n = 0:Nscale-1

    Mn = M_scale(n+1);

    for m = 0:Mn-1
        % Local time argument in Eq. (73)
        tLocal = q^n * (t - m/(q^n * W));

        % Evaluate g_tx at the scaled and shifted time points
        g_tx_nm = odssTransmitPulse(tLocal, q, W, T, pulseType);

        % ODSS subcarrier s_{m,n}(t), Eq. (73)
        s_mn = q^(n/2) * g_tx_nm;

        % Weighted superposition, Eq. (40)
        s_tx = s_tx + X_grid(n+1,m+1) * s_mn;
    end
end

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

    g_tx_n = odssTransmitPulse( ...
        tLocal, q, W, T, pulseType);

    s_n = q^(n/2) * g_tx_n;

    S_n = fft(s_n, Nfft);
    S_n = S_n(1:Nfft/2);

    S_n_dB = 20*log10( ...
        abs(S_n)/max(abs(S_n)) + eps);

    plot(f, S_n_dB, ...
        'DisplayName', sprintf('n = %d', n));
end

xlabel('Frequency (Hz)');
ylabel('Normalized magnitude (dB)');
title('Individual ODSS Scale Spectra');
xlim([0 B]);
ylim([-60 5]);
grid on;
legend('Location', 'best');