function [BER_LS, SNR_dB] = RX(folderName)
% Simulate the ODSS receiver over a pure AWGN channel.
% Input:
%   folderName - folder containing AWGN_Data.mat
% Outputs:
%   BER_LS - simulated BER corresponding to SNR_dB
%   SNR_dB - time-domain waveform SNR values
% Pure AWGN identity channel:
%   r(t) = s(t) + w(t)
% After canonical-dual matched filtering:
%   Yhat = X + W
% Therefore, no channel equalizer is required.

    %% 1. Load transmitter data
    dataFile = fullfile(folderName, 'TX_Data.mat');

    if exist(dataFile, 'file') ~= 2
        error('Cannot find TX_Data.mat in folder: %s', folderName);
    end

    load(dataFile, ...
        'rec_data', 'my_symbols', ...
        'tx_bits', 'x_data', 'X_data', 'signal_power', ...
        'num_paths', 'num_run', 'num_blocks', ...
        'QAM_order', 'bitsPerSymbol', ...
        'q', 'Nscale', 'B', 'W', 'T', 'Fs', 'Ts', ...
        'pulseType', ...
        'n_set', 'M_scale', 'M_tot', ...
        'pair_n', 'pair_m', ...
        'T_iMF', ...
        'transformRank', 'transformCondition', ...
        'Ns', 't', ...
        'G_tx', 'G_rx', 'Gram_tx', 'biorthErr', ...
        'chan_tau', 'chan_A', 'chan_a');

    %% 2. Noise covariance for theoretical QPSK BER
    % Unit-power time-domain noise:
    %   w ~ CN(0,I)
    % After delay-scale matched filtering:
    %   W = Ts*G_rx^H*w
    % Therefore:
    %   Cov(W) = Ts^2*G_rx^H*G_rx
    noiseCovUnitY = Ts^2*(G_rx'*G_rx);
    noiseCovUnitY = (noiseCovUnitY + noiseCovUnitY')/2;
    
    sourceNoiseCovUnit = T_iMF \ (noiseCovUnitY / T_iMF');
    sourceNoiseCovUnit = (sourceNoiseCovUnit + sourceNoiseCovUnit')/2;
    
    sourceNoiseVarianceUnit = real(diag(sourceNoiseCovUnit));
    
    if any(sourceNoiseVarianceUnit <= 0)
        error('Calculated source-domain noise variance is not positive.');
    end
    
    %% Average source-domain noise gain
    meanSourceNoiseGain = mean(sourceNoiseVarianceUnit);
    
    sourceNoiseSpread_dB = 10*log10(max(sourceNoiseVarianceUnit)/min(sourceNoiseVarianceUnit));
    
    fprintf('Mean source-domain noise gain = %.6e\n', meanSourceNoiseGain);
    
    fprintf('Source-domain noise variance min/max = %.6e / %.6e\n', ...
        min(sourceNoiseVarianceUnit), ...
        max(sourceNoiseVarianceUnit));
    
    fprintf('Source-domain noise spread = %.2f dB\n', sourceNoiseSpread_dB);

    %% Complete transmission matrix
    A_tx = G_tx*T_iMF;
    
    %% Complete receiver matrix
    B_rx = T_iMF \ (Ts*G_rx');
    
    %% Condition 1
    leftInverseError = norm( ...
        B_rx*A_tx-eye(M_tot), ...
        'fro')/sqrt(M_tot);

    %% Condition 2, calculated in two equivalent ways
    Cunit_from_B = B_rx*B_rx';
    
    % Avoid inv(); this calculates inv(A_tx'*A_tx)
    Cunit_from_A = ...
        (A_tx'*A_tx) \ eye(M_tot);
    
    covarianceIdentityError = norm( ...
        Cunit_from_B-Cunit_from_A, ...
        'fro') / norm(Cunit_from_A,'fro');
    
    Cunit = (Cunit_from_A+Cunit_from_A')/2;
    
    noiseVarianceUnit = real(diag(Cunit));
    meanNoiseGain = mean(noiseVarianceUnit);
    
    noiseGainSpread_dB = 10*log10( ...
        max(noiseVarianceUnit) / ...
        min(noiseVarianceUnit));
    
    offDiagonal = Cunit-diag(diag(Cunit));
    
    noiseCorrelationRatio = ...
        norm(offDiagonal,'fro') / norm(Cunit,'fro');
    
    iidNoiseError = norm( ...
        Cunit-meanNoiseGain*eye(M_tot), ...
        'fro') / ...
        norm(meanNoiseGain*eye(M_tot),'fro');
    
    fprintf('\nComplete operator diagnostics\n');
    fprintf('BA-I error                    = %.3e\n', ...
        leftInverseError);
    fprintf('BB^H versus inv(A^H A) error = %.3e\n', ...
        covarianceIdentityError);
    fprintf('Mean source noise gain        = %.6e\n', ...
        meanNoiseGain);
    fprintf('Noise-gain spread             = %.3f dB\n', ...
        noiseGainSpread_dB);
    fprintf('Noise correlation ratio       = %.3e\n', ...
        noiseCorrelationRatio);
    fprintf('IID-noise error               = %.3e\n', ...
        iidNoiseError);

    %% 6. SNR range
    % SNR is defined at the recovered source-QPSK domain:
    %
    %   SNR = Es_source/average source-domain noise variance
    %
    % qammod(...,'UnitAveragePower',true), so:
    %
    %   Es_source = 1
    %
    % Since:
    %
    %   average source noise variance
    %       = noisePower*meanSourceNoiseGain
    %
    % choose:
    %
    %   noisePower = 1/(SNR*meanSourceNoiseGain)
    SNR_dB = 0:3:20;

    numFrames = num_run*num_blocks;

    BER_LS = zeros(1,length(SNR_dB));
    SER_LS = zeros(1,length(SNR_dB));

    BER_theory = zeros(1,length(SNR_dB));
    BER_QPSK = zeros(1,length(SNR_dB));
    measuredSNR_LS = zeros(1,length(SNR_dB));

    %% 7. Loop over SNR
    for ss = 1:length(SNR_dB)

        BER_run = zeros(1,numFrames);
        SER_run = zeros(1,numFrames);

        measuredSNR_run = zeros(1,numFrames);

        %% 8. Loop over runs and blocks
        for frameIndex = 1:numFrames
            runIndex = ceil(frameIndex/num_blocks);
            blockIndex = frameIndex - (runIndex-1)*num_blocks;

            %% 8.1 Load transmitter data
            s_tx = rec_data(:,blockIndex,runIndex);

            x_true = x_data(:,blockIndex,runIndex);
            X_true = X_data(:,blockIndex,runIndex);

            bits_true = double(tx_bits(:,blockIndex,runIndex));
            symbols_true = double(my_symbols(:,blockIndex,runIndex));

            %% 8.2 Calculate Eq. (65) for the current channel realization
            H = calc_ODSS_H_matrix(chan_tau(:,runIndex), chan_A(:,runIndex), chan_a(:,runIndex), q, W, T, pulseType, t, M_scale);
            D = diag(diag(H));
            D = eye(M_tot);
            
            
            %% 8.3 Add complex time-domain AWGN
            SNR_linear = 10^(SNR_dB(ss)/10);
            
            % Unit-average-power source QPSK: Es = 1
            Es_source = 1;
            
            % Make the average source-domain noise variance equal to 1/SNR
            noisePower = Es_source/(SNR_linear*meanSourceNoiseGain);
            
            w_time = sqrt(noisePower/2).*(randn(Ns,1) + 1j*randn(Ns,1));
            
            y = s_tx + w_time;
            
             %% 8.4 ODSS delay-scale matched-filter demodulation
            Yhat_vector = complex(zeros(M_tot,1));
            
            for idx = 1:M_tot
                template_nm = G_rx(:,idx);
                g_rx_nm_fun = @(u) interp1(t,template_nm,u,'linear',0);
            
                Yhat_vector(idx) = delay_scale_ambiguity(y,t,g_rx_nm_fun,0,1);
            end
            
            %% 8.5 One-tap MMSE equalization, Eq. (69)
            sigmaW2 = noisePower*real(trace(noiseCovUnitY))/M_tot;
            
            Zhat_vector = D'*((D*D' + sigmaW2*eye(M_tot))\Yhat_vector);
            
            %% 8.6 Inverse Mellin-Fourier transform, Eq. (70)
            x_hat = T_iMF\Zhat_vector;
            
            %% Actual decision-domain SNR for diagnostics
            mmseGain = 1/(1+sigmaW2);
            
            x_signal = mmseGain*x_true;
            x_noise_actual = x_hat-x_signal;
            
            measuredSNR_run(frameIndex) = 10*log10(mean(abs(x_signal).^2)/mean(abs(x_noise_actual).^2));

            %% 8.6 QAM detection
            symbols_hat = qamdemod(x_hat, QAM_order, 'gray', 'UnitAveragePower', true);
            symbols_hat = double(symbols_hat(:));
            bits_hat_matrix = de2bi(symbols_hat, bitsPerSymbol, 'left-msb');
            bits_hat = reshape(bits_hat_matrix.',[],1);

            %% 8.7 BER
            numBitErrors = nnz(bits_hat ~= bits_true);
            BER_run(frameIndex) = numBitErrors/numel(bits_true);

            %% 8.8 SER
            numSymbolErrors = nnz(symbols_hat ~= symbols_true);
            SER_run(frameIndex) = numSymbolErrors/M_tot;


            %% 8.10 Save final-frame variables
            if ss == length(SNR_dB) && frameIndex == numFrames
                final_y = y;
                final_w_time = w_time;

                final_Yhat_vector = Yhat_vector;
                final_Zhat_vector = Zhat_vector;

                final_X_true = X_true;
                final_x_soft = x_hat;
                final_x_true = x_true;

                final_symbols_hat = symbols_hat;
                final_bits_hat = bits_hat;
            end
        end


        %% 9. Average results
        BER_LS(ss) = mean(BER_run);
        SER_LS(ss) = mean(SER_run);
        measuredSNR_LS(ss) = mean(measuredSNR_run);
        
        if QAM_order == 4
            % Current ODSS analytical BER.
            % Each recovered source symbol may have a different noise variance.
            sourceNoiseVariance = noisePower*sourceNoiseVarianceUnit;
        
            BER_theory(ss) = mean(qfunc(sqrt(1./sourceNoiseVariance)));
        
            % Standard Gray-QPSK BER under equal white-noise variance
            BER_QPSK(ss) = qfunc(sqrt(SNR_linear));
        else
            BER_theory(ss) = NaN;
            BER_QPSK(ss) = NaN;
        end
        
        fprintf(['SNR = %5.1f dB, measured decision SNR = %6.2f dB, ', 'BER = %.6e\n'], ...
            SNR_dB(ss),measuredSNR_LS(ss),BER_LS(ss));

    end

    %% 10. Plot BER curve
    figure;
    
    if QAM_order == 4
        semilogy(SNR_dB,BER_LS,'b-o', ...
                 SNR_dB,BER_theory,'m+--', ...
                 SNR_dB,BER_QPSK,'r*--', ...
                 'LineWidth',1.5);
    
        legend('ODSS simulation BER', ...
               'ODSS analytical BER', ...
               'Standard QPSK BER', ...
               'Location','southwest');
    else
        semilogy(SNR_dB,BER_LS,'b-o','LineWidth',1.5);
        legend('ODSS simulation BER','Location','southwest');
    end
    
    grid on;
    ylim([1e-6 1]);
    
    xlabel('Decision-domain E_s/N_0 (dB)');
    ylabel('Bit Error Rate');
    
    title(sprintf('ODSS over pure AWGN, q = %.2f, N_{scale} = %d', q,Nscale));

    %% 11. Plot final constellation
    figure;
    plot(real(final_x_soft), imag(final_x_soft), 'o');
    hold on;
    plot(real(final_x_true), imag(final_x_true), 'kx', 'LineWidth', 1.5);
    grid on;
    axis equal;
    xlabel('In-phase');
    ylabel('Quadrature');
    title(sprintf( 'ODSS recovered symbols, SNR = %.1f dB', SNR_dB(end)));
    legend( 'Recovered symbols', 'Transmitted symbols', 'Location', 'best');

    %% 12. Save receiver results
    rxFormatVersion = 1;
    rxTag = sprintf('RX_q%g_Nscale%d', q, Nscale);
    rxDataFile = fullfile(folderName, [rxTag '.mat']);
    save(rxDataFile, ...
        'rxFormatVersion', 'SNR_dB', 'BER_LS', 'SER_LS', 'BER_theory', 'measuredSNR_LS', ...
        'q', 'Nscale', 'M_tot', 'QAM_order', 'num_run', 'num_blocks', ...
        'numFrames', 'noiseCovUnitY', 'sourceNoiseCovUnit', 'sourceNoiseVarianceUnit', ...
        'biorthErr', 'final_y', 'final_w_time', ...
        'final_Yhat_vector', 'final_Zhat_vector', 'final_X_true', ...
        'final_x_soft', 'final_x_true', 'final_symbols_hat', 'final_bits_hat', 'BER_QPSK', ...
        'meanSourceNoiseGain', 'sourceNoiseSpread_dB');

    fprintf('\nReceiver data saved successfully:\n%s\n', rxDataFile);
end