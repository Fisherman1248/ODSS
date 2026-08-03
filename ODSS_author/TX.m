function [folderName] = TX(q)
% Generate transmitter data for the ODSS pure-AWGN simulation.
%
% Input:
%   q          - ODSS scale factor
%
% Output:
%   folderName - folder containing AWGN_Data.mat
%
% Notes:
%   1. No CP or zero padding is used.
%   2. No multipath delay or Doppler is applied.
%   3. No AWGN is added in this function.
%      The receiver adds noise according to the required SNR.
%   4. rec_data contains the noiseless baseband ODSS waveform.

    %% 1. ODSS parameters
    Nscale = 6;
    B = 1280;
    T = 1.9;
    Fs = 10240;
    Ts = 1/Fs;

    pulseType = 'phydyas';
    QAM_order = 4;

    if q <= 1
        error('q must be greater than 1.');
    end

    W = B*(q-1)/(q^Nscale-1);

    %% 2. Monte Carlo parameters
    num_paths = 1;
    num_run = 100;
    num_blocks = 1;

    bitsPerSymbol = log2(QAM_order);

    if mod(bitsPerSymbol,1) ~= 0
        error('QAM_order must be a power of two.');
    end

    %% 3. Valid ODSS index mapping
    % pair_n(idx) and pair_m(idx) contain the mathematical
    % ODSS index (n,m) associated with vector index idx.
    n_set = 0:Nscale-1;
    M_scale = floor(q.^n_set);
    M_tot = sum(M_scale);

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

    %% 4. Construct the discrete inverse Mellin-Fourier matrix
    % X_vector = T_iMF*x_vector
    %
    % This implements Eq. (38):
    %
    % X[n,m] = q^(-n/2)/Nscale *
    %          sum_k sum_l x[k,l]/M(k) *
    %          exp(j*2*pi*(m*l/M(k) - n*k/Nscale))
    T_iMF = complex(zeros(M_tot,M_tot));

    for outIdx = 1:M_tot
        n = pair_n(outIdx);
        m = pair_m(outIdx);

        for inIdx = 1:M_tot
            k = pair_n(inIdx);
            l = pair_m(inIdx);

            Mk = M_scale(k+1);

            T_iMF(outIdx,inIdx) = q^(-n/2)/(Nscale*Mk) * exp(1j*2*pi*(m*l/Mk - n*k/Nscale));
        end
    end

    transformRank = rank(T_iMF);
    transformCondition = cond(T_iMF);

    %% 5. Construct the sampled ODSS transmit basis
    Ns = round(T*Fs);

    % Useful waveform interval: [0,T)
    t = (0:Ns-1).'/Fs;

    G_tx = complex(zeros(Ns,M_tot));

    for col = 1:M_tot
        n = pair_n(col);
        m = pair_m(col);

        % Argument of the prototype pulse:
        % q^n(t - m/(q^n W))
        tLocal = q^n*(t - m/(q^n*W));

        G_tx(:,col) = q^(n/2) * odssTransmitPulse(tLocal,q,W,T,pulseType);
    end
    

    %% 6. Construct the canonical dual receive basis
    Gram_tx = Ts*(G_tx'*G_tx);

    if rcond(Gram_tx) < 1e-12
        warning(['The transmit Gram matrix is poorly conditioned. ', ...
            'rcond(Gram_tx) = %.3e'],rcond(Gram_tx));
    end

    G_rx = G_tx/Gram_tx;

    biorthErr = norm(Ts*(G_rx'*G_tx)-eye(M_tot),'fro')/sqrt(M_tot); % Ts * G_rx^H * G_tx = I

    %% 7. Preallocate Monte Carlo data
    % rec_data(:,ii,zz):
    %   received noiseless waveform for block ii and run zz
    %
    % my_symbols(:,ii,zz):
    %   integer QAM symbol labels
    %
    % tx_bits(:,ii,zz):
    %   source bits
    %
    % x_data(:,ii,zz):
    %   source-domain QAM vector x[k,l]
    %
    % X_data(:,ii,zz):
    %   inverse Mellin-Fourier output X[n,m]

    rec_data = complex(zeros(Ns,num_blocks,num_run));

    my_symbols = zeros(M_tot,num_blocks,num_run,'uint16');

    tx_bits = zeros(M_tot*bitsPerSymbol,num_blocks,num_run,'uint8');

    x_data = complex(zeros(M_tot,num_blocks,num_run));

    X_data = complex(zeros(M_tot,num_blocks,num_run));

    signal_power = zeros(num_blocks,num_run);

    %% 8. Pure-AWGN channel parameters
    % These variables are retained to keep the saved-data structure
    % consistent with the previous OHFM implementation.
    %
    % Pure AWGN identity channel:
    %   delay             = 0
    %   channel gain      = 1
    %   Doppler parameter = 0
    chan_tau = zeros(num_paths,num_run);
    chan_A = ones(num_paths,num_run);
    chan_a = zeros(num_paths,num_run);

    %% 9. Transmitter Monte Carlo loop
    for zz = 1:num_run

        for ii = 1:num_blocks
            bits = randi([0 1], M_tot*bitsPerSymbol, 1); % Generate source bits

            % Convert bits to integer QAM symbols
            symbols = bi2de(reshape(bits,bitsPerSymbol,[]).', 'left-msb');

            x_vector = qammod( symbols, QAM_order, 'gray', 'UnitAveragePower',true); % QAM modulation
            X_vector = T_iMF*x_vector; % Discrete inverse Mellin-Fourier transform

            s_tx = G_tx*X_vector; % ODSS waveform modulation

            % Pure AWGN identity channel
            % Noise is not added here.
            y = s_tx;

            % Store generated data
            tx_bits(:,ii,zz) = uint8(bits);
            my_symbols(:,ii,zz) = uint16(symbols);

            x_data(:,ii,zz) = x_vector;
            X_data(:,ii,zz) = X_vector;

            rec_data(:,ii,zz) = y;

            signal_power(ii,zz) = mean(abs(s_tx).^2);
        end

        if mod(zz,10) == 0 || zz == 1 || zz == num_run
            fprintf('Generated run %d/%d\n',zz,num_run);
        end
    end

    %% 10. Generate output folder
    folderName = sprintf('ODSS_AWGN_q%g_Nscale%d', q,Nscale);

    if ~exist(folderName,'dir')
        mkdir(folderName);
    end

    %% 11. Save transmitter data
    txFormatVersion = 3;

    dataFile = fullfile(folderName,'TX_Data.mat');

    save(dataFile, ...
        'txFormatVersion', 'rec_data', 'my_symbols', ...
        'tx_bits', 'x_data', 'X_data', 'signal_power', ...
        'num_paths', 'num_run', 'num_blocks', 'QAM_order', 'bitsPerSymbol', ...
        'q', 'Nscale', 'B', 'W', 'T', 'Fs', 'Ts', 'pulseType', ...
        'n_set', 'M_scale', 'M_tot', 'pair_n', 'pair_m', ...
        'T_iMF', 'transformRank', 'transformCondition', 'Ns', 't', ...
        'G_tx', 'G_rx', 'Gram_tx', 'biorthErr', ...
        'chan_tau', 'chan_A', 'chan_a', '-v7.3');

    %% 12. Display transmitter diagnostics
    fprintf('\n');
    fprintf('ODSS transmitter data saved successfully.\n');
    fprintf('File: %s\n',dataFile);
    fprintf('q = %.4f\n',q);
    fprintf('Nscale = %d\n',Nscale);
    fprintf('M_tot = %d\n',M_tot);
    fprintf('Waveform samples per block = %d\n',Ns);
    fprintf('Number of runs = %d\n',num_run);
    fprintf('Number of blocks per run = %d\n',num_blocks);
    fprintf('rank(T_iMF) = %d/%d\n',transformRank,M_tot);
    fprintf('cond(T_iMF) = %.3e\n',transformCondition);
    fprintf('Biorthogonality error = %.3e\n',biorthErr);
    fprintf('Average waveform power = %.6e\n', mean(signal_power(:)));

end