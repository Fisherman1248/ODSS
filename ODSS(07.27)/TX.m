%% ODSS transmitter
% Generate an ODSS waveform and save everything needed by RX.m.
%
% Output:
%   TX_data.mat (saved next to this script)
%
% MATLAB R2016b or newer is required for local functions in scripts.

clear; clc; close all;

%% 1. Transmitter parameters
q = 1.5;
Nscale = 6;
B = 1280;
W = B*(q-1)/(q^Nscale-1);
T = 1.9;
Fs = 10240;
Ts = 1/Fs;
pulseType = 'phydyas';
QAM_order = 4;
SNRdB_curve = 0:2:40;

n_set = 0:Nscale-1;
M_scale = floor(q.^n_set);
M_tot = sum(M_scale);

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

%% 3. Generate the source bits and QAM symbols
bitsPerSymbol = log2(QAM_order);
tx_bits = randi([0,1], M_tot*bitsPerSymbol, 1);
x_vector = qammod(tx_bits, QAM_order, 'InputType', 'bit', 'UnitAveragePower', true);

%% 4. Map the serial QAM symbols to x[k,l]
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

%% 6. Convert X_grid to X_vector
X_vector = complex(zeros(M_tot,1));
symbolPointer = 1;
for n = 0:Nscale-1
    Mn = M_scale(n+1);
    symbolIndices = symbolPointer:symbolPointer+Mn-1;
    X_vector(symbolIndices) = X_grid(n+1,1:Mn).';
    symbolPointer = symbolPointer + Mn;
end

%% 7. Construct T_iMF for receiver-side inverse transformation
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

%% 8. Build the sampled ODSS synthesis matrix
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

%% 9. ODSS waveform modulation
s_tx = G_tx*X_vector;

%% 10. Save transmitter output and receiver metadata
% G_tx is intentionally not stored because RX.m reconstructs its own
% receive basis from the saved waveform parameters.
txFormatVersion = 1;
txDataFile = fullfile(fileparts(mfilename('fullpath')), 'TX_data.mat');

save(txDataFile, ...
    'txFormatVersion', ...
    'q', 'Nscale', 'B', 'W', 'T', 'Fs', 'Ts', 'SNRdB_curve', ...
    'pulseType', 'QAM_order', 'bitsPerSymbol', ...
    'n_set', 'M_scale', 'M_tot', 'Ns', 't', ...
    'pair_n', 'pair_m', 'T_iMF', ...
    'tx_bits', 'x_vector', 'x_grid', ...
    'X_vector', 'X_grid', 's_tx', 'G_tx', 'G_rx', ...
    'loopMatrixError');

fprintf('Transmitter data saved successfully:\n%s\n', txDataFile);
fprintf('Saved waveform samples: %d\n', numel(s_tx));
fprintf('Saved source bits: %d\n', numel(tx_bits));
