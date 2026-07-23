%% ODSS AWGN smoke test -- corrected from your original script
% Requires your existing odssTransmitPulse.m on the MATLAB path.
% This version intentionally removes delay/scale channel and Eq. (65)
% one-tap equalization. It first verifies the sampled AWGN closed loop.

clc; clear; close all;
rng(1);

%% 1. Simulation setup
q = 2;
Nscale = 7;
B = 1280;
W = B*(q-1)/(q^Nscale-1);
T = 1.9;
oversampling = 8;
pulseType = 'phydyas_2';
QAM_order = 4;
SNRdB = 20;

Fs = oversampling * B;
Ts = 1/Fs;
t = (0:Ts:T-Ts).';                 % column vector

n_set = 0:Nscale-1;
M_scale = floor(q.^n_set);
M_tot = sum(M_scale);

fprintf('q = %.1f, Nscale = %d, M_tot = %d\n', q, Nscale, M_tot);
fprintf('W = %.6f Hz, T = %.6f s, Fs = %.1f Hz\n', W, T, Fs);

%% 2. Generate QAM symbols
bitsPerSymbol = log2(QAM_order);
tx_bits = randi([0 1], M_tot*bitsPerSymbol, 1);
x_vector = qammod(tx_bits, QAM_order, ...
    'InputType','bit','UnitAveragePower',true);

%% 3. Map x_vector to x[k,l]
x_grid = complex(zeros(Nscale, M_scale(end)));
symbolPointer = 1;
for k = 0:Nscale-1
    Mk = M_scale(k+1);
    idx = symbolPointer:symbolPointer+Mk-1;
    x_grid(k+1,1:Mk) = x_vector(idx).';
    symbolPointer = symbolPointer + Mk;
end

%% 4. Forward inverse-Mellin/Fourier ODSS transform
% This is exactly the transform implemented by your original Eq. (38) loop.
X_grid = complex(zeros(Nscale, M_scale(end)));
for n = 0:Nscale-1
    Mn = M_scale(n+1);
    for m = 0:Mn-1
        temp = 0;
        for k = 0:Nscale-1
            Mk = M_scale(k+1);
            l = 0:Mk-1;
            phase = exp(1j*2*pi*(m*l/Mk - n*k/Nscale));
            temp = temp + sum(x_grid(k+1,1:Mk).*phase)/Mk;
        end
        X_grid(n+1,m+1) = q^(-n/2)*temp/Nscale;
    end
end

% Stack X[n,m] using the same n,m ordering as the waveform bank.
X_vector = complex(zeros(M_tot,1));
symbolPointer = 1;
for n = 0:Nscale-1
    Mn = M_scale(n+1);
    idx = symbolPointer:symbolPointer+Mn-1;
    X_vector(idx) = X_grid(n+1,1:Mn).';
    symbolPointer = symbolPointer + Mn;
end

%% 5. Construct the exact matrix of the implemented forward transform
T_iMF = complex(zeros(M_tot,M_tot));
outIndex = 1;
for n = 0:Nscale-1
    Mn = M_scale(n+1);
    for m = 0:Mn-1
        inIndex = 1;
        for k = 0:Nscale-1
            Mk = M_scale(k+1);
            for l = 0:Mk-1
                T_iMF(outIndex,inIndex) = ...
                    q^(-n/2)/(Nscale*Mk) * ...
                    exp(1j*2*pi*(m*l/Mk - n*k/Nscale));
                inIndex = inIndex + 1;
            end
        end
        outIndex = outIndex + 1;
    end
end

forwardMatrixError = norm(X_vector-T_iMF*x_vector)/max(norm(X_vector),eps);
transformOnlyError = norm((T_iMF\X_vector)-x_vector)/max(norm(x_vector),eps);

fprintf('\nTransform diagnostics\n');
fprintf('rank(T_iMF) = %d / %d\n',rank(T_iMF),M_tot);
fprintf('cond(T_iMF) = %.3e\n',cond(T_iMF));
fprintf('loop-versus-matrix error = %.3e\n',forwardMatrixError);
fprintf('transform-only inverse error = %.3e\n',transformOnlyError);

%% 6. Build sampled ODSS transmit waveform bank G_tx
% Column i is one sampled waveform associated with one X[n,m].
G_tx = complex(zeros(length(t),M_tot));
symbolPointer = 1;
for n = 0:Nscale-1
    Mn = M_scale(n+1);
    for m = 0:Mn-1
        tLocal = q^n*(t - m/(q^n*W));
        g_tx_nm = odssTransmitPulse(tLocal,q,W,T,pulseType);
        s_mn = q^(n/2)*g_tx_nm;
        G_tx(:,symbolPointer) = s_mn(:);
        symbolPointer = symbolPointer + 1;
    end
end

s_tx = G_tx*X_vector;

%% 7. Construct a sampled biorthogonal receive bank G_rx
% Ts*G_rx^H*G_tx = I.
Gram = Ts*(G_tx'*G_tx);
if rcond(Gram) < 1e-13
    error('G_tx Gram matrix is numerically singular.');
end
G_rx = G_tx/Gram;

biorthError = norm(Ts*(G_rx'*G_tx)-eye(M_tot),'fro')/sqrt(M_tot);

fprintf('\nWaveform diagnostics\n');
fprintf('samples = %d\n',length(s_tx));
fprintf('rank(G_tx) = %d / %d\n',rank(G_tx),M_tot);
fprintf('cond(Ts*G_tx^H*G_tx) = %.3e\n',cond(Gram));
fprintf('biorthogonality error = %.3e\n',biorthError);

%% 8. Noiseless waveform closed-loop test
Y0_vector = Ts*(G_rx'*s_tx);
x0_hat = T_iMF\Y0_vector;

Y0Error = norm(Y0_vector-X_vector)/max(norm(X_vector),eps);
x0Error = norm(x0_hat-x_vector)/max(norm(x_vector),eps);

fprintf('\nNoiseless closed-loop\n');
fprintf('||Y0-X||/||X|| = %.3e\n',Y0Error);
fprintf('||x0-x||/||x|| = %.3e\n',x0Error);

%% 9. Pure AWGN channel
% IMPORTANT: no random h, no interpolation, no tau, and no alpha.
rs = s_tx;
Ps = mean(abs(rs).^2);
Pn = Ps*10^(-SNRdB/10);
w = sqrt(Pn/2)*(randn(size(rs))+1j*randn(size(rs)));
r = rs+w;

measuredSNR = 10*log10(mean(abs(rs).^2)/mean(abs(w).^2));

fprintf('\nPure AWGN channel\n');
fprintf('signal power = %.6e\n',Ps);
fprintf('noise power = %.6e\n',Pn);
fprintf('target SNR = %.2f dB\n',SNRdB);
fprintf('measured SNR = %.2f dB\n',measuredSNR);

%% 10. ODSS receiver
A = G_tx * T_iMF;

% w 中每个复采样点的噪声方差
sigmaW2 = Pn;

x_hat_vector = ...
    (A' * A + sigmaW2 * eye(M_tot)) \ ...
    (A' * r);

XError = norm(Zhat_vector-X_vector)/max(norm(X_vector),eps);
dataSymbolError = norm(x_hat_vector-x_vector)/max(norm(x_vector),eps);

fprintf('\nReceiver diagnostics\n');
fprintf('Recovered X relative error = %.3e\n',XError);
fprintf('Recovered data relative error = %.3e\n',dataSymbolError);

%% 11. QAM decision
rx_bits = qamdemod(x_hat_vector,QAM_order, ...
    'OutputType','bit','UnitAveragePower',true);
rx_bits = rx_bits(:);

numBitErrors = sum(rx_bits~=tx_bits);
BER = numBitErrors/numel(tx_bits);

fprintf('\nReceiver results\n');
fprintf('Bit errors: %d / %d\n',numBitErrors,numel(tx_bits));
fprintf('BER: %.6e\n',BER);

%% 12. Constellations
figure;
subplot(1,2,1);
scatter(real(x_vector),imag(x_vector),25,'filled');
grid on; axis equal;
xlabel('In-phase'); ylabel('Quadrature');
title('Transmitted QAM symbols');

subplot(1,2,2);
scatter(real(x_hat_vector),imag(x_hat_vector),25,'filled');
grid on; axis equal;
xlabel('In-phase'); ylabel('Quadrature');
title(sprintf('Recovered QAM symbols, SNR = %.1f dB',SNRdB));
