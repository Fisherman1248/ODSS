% ------------ Spectral efficiency evaluation: CP-OFDM, VBMC --------------
%
% This code computes the maximum spectral efficiency of various
% communication schemes (CP-OFDM, VBMC) over a delay-scale spread wideband
% channel.
%
% A delay and scale spread channel with multiple paths is generated.
% Correlation rake-receiver is implemented as part of the receiver
% processing. Perfect channel state information is assumed at the receiver
% side.
%
%Authors : (1) Arunkumar K. P.
%          (2) Chandra R. Murthy
%          (3) Dr. P. Muralikrishna
%Address : (1) Ph.D. Scholar,
%              Signal Processing for Communications Lab, ECE Department,
%              Indian Institute of Science, Bangalore, India-560 012.
%          (2) Professor,
%              Electrical Communication Engineering (ECE) Department,
%              Indian Institute of Science, Bangalore, India-560 012.
%          (3) Scientist,
%              Naval Physical Oceanographic Laboratory, Kochi - 682021
%Email   : arunkumar@iisc.ac.in
%
%Revision History
% Version : 1.0
% Created : 04-01-2022
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% This is a script under development.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

defpaper
clear;
clc;
close all;

%% Channel parameters
tauMax = 10e-3; %max delay
alphaMax = 1;%1.01;%1.01; %max time-scale
Np = 10; %number of paths from Tx to Rx
SNRvalues = inf;%0:3:30;%SNR values, in dB

%% System parameters
fL = 10e3;  %start frequency, in Hz, of the communication frequency band
fH = 20e3;  %end frequency, in Hz, of the communication frequency band
W = fH - fL;%bandwidth, in Hz, of the communication frequency band

%% CP-OFDM system parameters
N = 128;   %number of symbols
Fs = W;     %sampling frequency (used in CP-OFDM system)
Ts = 1/Fs;  %sampling interval
T = N*Ts;   %symbol duration (CP-OFDM)
L = ceil(tauMax*Fs); %CP length

%% Monte Carlo Trials
Ntrials = 100;% # of Monte Carlo trials to estimate BER
if (length(SNRvalues)==1) && eq(SNRvalues,inf)
    Ntrials = 1; %for infinite SNR need to loop only once
end

NBITERRMIN = 100; % # of bit errors to wait for, before closing the Monte Carlo run
NtrialMin = 10;

bitErrFullCSIR = zeros(size(SNRvalues)); %BER of the decoder using Full CSI
bitErrPartCSIR = zeros(size(SNRvalues)); %BER of the decoder using partial, imperfect CSI (only diagonal of Hmat_, ignoring time scale)

CH = zeros(Ntrials, length(SNRvalues)); %maximum spectral efficiency

PrSignal = zeros(Ntrials, length(SNRvalues));
PrNoise = zeros(Ntrials, length(SNRvalues));

PrSignal_ = zeros(Ntrials, length(SNRvalues));
PrNoise_ = zeros(Ntrials, length(SNRvalues));

PyNoise = zeros(Ntrials, length(SNRvalues));
PySignal = zeros(Ntrials, length(SNRvalues));

displayRate = 30/60; %results are displayed at displayRate (minutes)
tDispl  = tic; % start the timer for displaying results

for iSNR = 1:length(SNRvalues)
    
    SNR = SNRvalues(iSNR) - 10*log10(Fs/W);
    
    for trial = 1:Ntrials
        
        disp(['SNR = ' num2str(SNR) ' dB, trial = ' num2str(trial)...
            ' of ' num2str(Ntrials) ' trials.']);
        
        %% Modulator (mounts data symbols)
        x = (randi(2,[1 N])-1-1/2)*2;%(randi(2,[1 Nb])-3/2)*2; %bits to be transmitted
        s = OFDM_modulation(1, N, L, x);%transmitted signal -- OFDM modulator output
        
        %         figure(6), clf,
        %         subplot(211), stem(x), grid on,
        %         xlabel('symbol index'), ylabel('symbol value'), title('transmitted symbols')
        %         subplot(212), plot(t, real(s)), grid on, hold on, plot(t, imag(s))
        %         xlabel('t (s)'), ylabel('s(t)'), title('transmitted signal')
        
        %% Delay-scale spread channel simulation
        taup = [0; rand(Np-1,1)*tauMax]; %arrival delays for each path -- perfect synchronisation is assumed, so the first path delay can be set to 0        
        taup = sort(taup,'ascend');%arrange delays in the order of their arrivals
        
        lAlpha = alphaMax - 1/alphaMax;%max. scale spread
        alphap = lAlpha*rand(Np,1) + 1/alphaMax; %time-scales for each path
        
        hp = (randn(Np,1) + 1i*randn(Np,1))/sqrt(2);%channel tap coefficients %[1; 2*(rand(Np-1,1)-0.5)];
        
        %         figure(8), clf, %plots channel parameters
        %         subplot(311), stem(taup), grid on,
        %         xlabel('path index, p'), ylabel('\tau(p)'),
        %         title('channel path delays')
        %         subplot(312), stem(alphap), grid on,
        %         xlabel('path index, p'), ylabel('\alpha(p)'),
        %         title('channel path time-scales')
        %         subplot(313), stem(taup, real(hp)), grid on, hold on, stem(taup, imag(hp)),
        %         xlabel('path delay, \tau(p)'), ylabel('h(p)'),
        %         title('channel path amplitudes'), legend('Re', 'Im')
        
        %% CSI
        % %% (1) Full CSI (NxN channel matrix, considering the effect of both delay and scale)
        % %% (2) Partial CSI (ignoring the effect of time-scale)
        H = zeros(N, N);% Full CSI: NbxNb channel matrix
        for iSymbol = 1:N
            bp = zeros(1, N);%all-zero transmission execpt for a pilot bit placed at the nb'th bit location
            bp(iSymbol) = 1; %pilot bit at the nb'th location - for perfect CSI extraction
            sp = OFDM_modulation(1, N, L, bp); %modulator output for pilot transmission
            
            rp = channel( sp, Fs, taup, alphap, hp, inf, 'trim', 100*Fs ); %replace alphap with 1's to avoid resampling
            %             rp  = channel( sp, Fs, taup, ones(size(alphap)), hp, inf ); %replace alphap with 1's to avoid resampling
            
            yp = OFDM_demodulation( 1, N, L, rp ); % OFDM demodulation
            
            H(:, iSymbol) = yp(:);
            
            % Display progress
            % fprintf(' %d', iSymbol);
            disp(['iSymbol = ' num2str(iSymbol)]);
        end
        %         H = diag(Hdiag);
        
        ICImax = max(max(20*log10(...
            bsxfun(@times, abs(H - diag(diag(H))), 1./diag(abs(H))')...
            )));
        
        disp(['Max. ICI (dB): ' num2str(ICImax)]);
        
        CH(trial, iSNR) = log2(abs(det(eye(N) + H*H')))/N; %spectral efficiency at 0 dB
        
        
        %% Receiver
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % 1) Front end signal received
        [ r, Ps, Pn, rs, w ] = channel( s, Fs, taup, alphap, hp, SNR, 'trim', 100*Fs );
        PrSignal(trial, iSNR) = Ps;    PrNoise(trial, iSNR) = Pn;
        
        %         [ r_, Ps_, Pn_] = channel( s, Fs, taup, ones(size(alphap)), hp, SNR );
        %         PrSignal_(trial, iSNR) = Ps_;    PrNoise_(trial, iSNR) = Pn_;
        
        %         figure(9), clf,%plots the transmitted signal and the received signal
        %         subplot(211), plot(t, real(s)), grid on, hold on, plot(t, imag(s)),
        %         xlabel('t (s)'), ylabel('s(t)'), legend('Re', 'Im')
        %         title('Transmitted signal')
        %         tr = (0:length(r)-1)/Fs; %sampling times for the received signal
        %         subplot(212), plot(tr, real(r)), grid on, hold on, plot(tr, imag(r))
        %         xlabel('t (s)'), ylabel('r(t)'), legend('Re', 'Im')
        %         title('Received signal after passing s(t) through delay-scale channel')
        
        % 2) Receiver side processing: OFDM demodulation
        y = OFDM_demodulation(1,N,L,r);%,diag(ones(M,1)));
        
        % 3) Receiver side processing
        % a) decode using full CSI
        Heff = H;
        x_est = (Heff'*Heff + Pn*eye(N))\(Heff'*y(:)); %MMSE decoder
        %         x_est = (Heff*Tp)\y(:); %LS decoder
        x_est = reshape(x_est, 1, N);
        
        figure(12),
        subplot(211), cla, stem(x), hold on, stem(real(x_est),'*'), grid on, stem(imag(x_est),'.')
        xlabel('bit index'), ylabel('bit value'), ylim([-2 2])
        legend('Transmitted bit', 'Received (soft) bit: Real', 'Received (soft) bit: Imag')
        title('Symbol recovery using full CSI')
        
        bitErrFullCSIR(iSNR) = ((trial-1)*bitErrFullCSIR(iSNR) +...
            sum(ne(x>0, real(x_est)>0))/length(x) )/trial; %symbol slicing and bit error evaluation
        
        disp(['BER (full CSIR): ' num2str(bitErrFullCSIR(iSNR))]);
        
        % b) decode using single tap equalizer
        Hd = diag(diag(H));
        xs_est = (Hd'*Hd + Pn*eye(N))\(Hd'*y(:)); %1-tap MMSE decoder
        xs_est = reshape(xs_est, 1, N);
        
        figure(12),
        subplot(212), cla, stem(x), hold on, stem(real(xs_est),'*'), grid on, stem(imag(xs_est),'.')
        xlabel('bit index'), ylabel('bit value'), ylim([-2 2])
        legend('Transmitted bit', 'Received (soft) bit: Real', 'Received (soft) bit: Imag')
        title('Symbol recovery using single tap equalizer')
        pause(0.01);
        
        bitErrPartCSIR(iSNR) = ((trial-1)*bitErrPartCSIR(iSNR) +...
            sum(ne(x>0, real(xs_est)>0))/length(x) )/trial; %symbol slicing and bit error evaluation
        
        disp(['BER (partial CSIR, single tap equalizer): ' num2str(bitErrPartCSIR(iSNR))]);
        
        %% Quit Monte-Carlo loop if atleast 100 bit errors are registered for each active decoder
        minBitErr =  min([bitErrFullCSIR(iSNR) bitErrPartCSIR(iSNR)]*trial*N);
        if (minBitErr > NBITERRMIN) && (trial > NtrialMin)
            disp(' ');
            disp([num2str(trial) ': ' num2str([bitErrFullCSIR(iSNR) bitErrPartCSIR(iSNR)])...
                ' (SNR = ' num2str(SNR) ' dB)' ' [bitErrs = ' num2str(minBitErr) ']']);
            disp(['Monte-Carlo loop terminated since total bit error count > '...
                num2str(NBITERRMIN) ' for all active decoders and # of trials='...
                num2str(trial) ' > ' num2str(NtrialMin)]);
            %Update the plot before proceeding to next SNR
            figure(18), cla
            semilogy(SNRvalues(1:iSNR), bitErrFullCSIR(1:iSNR), '-*'), hold on, grid on
            semilogy(SNRvalues(1:iSNR), bitErrPartCSIR(1:iSNR), '-o'), hold on, grid on
            xlabel('SNR (dB)'), ylabel('BER'),
            legend(['CP-OFDM (T = ' num2str(T) ', Full CSI)'], ...
                ['CP-OFDM (T = ' num2str(T) ', Partial CSI)'])
            title([num2str(trial) '/' num2str(Ntrials) ' trials, bit errors (best scheme) = ' num2str(minBitErr)])
            
            figure(99), cla
            plot(1:trial, PrSignal(1:trial,iSNR), '-*'), hold on, grid on
            plot(1:trial, PrSignal_(1:trial,iSNR), '-o'), hold on, grid on
            plot(1:trial, ones(trial,1)*Np, 'k'), %overlay expected signal power
            xlabel('trial #'), ylabel('P_{signal}'), legend('Full CSI', 'Partial CSI')
            title(['Mean Rcvd Sig Pow: [' num2str(mean(PrSignal(1:trial,iSNR)))...
                ', ' num2str(mean(PrSignal_(1:trial,iSNR))) ']'])
            pause(0.1) %waiting a bit to get the plot window updated
            break;
        end
        
        %% Display results (in the run time)
        if toc(tDispl) > displayRate*60
            %Save the workspace into a mat file every saveRate minutes
            tDispl  = tic; % restart the timer for saving results
            disp(' ');
            disp([num2str(trial) ': ' num2str([bitErrFullCSIR(iSNR) bitErrPartCSIR(iSNR)])...
                ' (SNR = ' num2str(SNR) ' dB)' ' [bitErrs = ' num2str(minBitErr) ']']);
            %Update the plot every displayRate minutes
            figure(18), cla
            semilogy(SNRvalues(1:iSNR), bitErrFullCSIR(1:iSNR), '-*'), hold on, grid on
            semilogy(SNRvalues(1:iSNR), bitErrPartCSIR(1:iSNR), '-o'), hold on, grid on
            xlabel('SNR (dB)'), ylabel('BER'),
            legend(['CP-OFDM (T = ' num2str(T) ', Full CSI)'],...
                ['CP-OFDM (T = ' num2str(T) ', Partial CSI)'])
            title([num2str(trial) '/' num2str(Ntrials) ' trials, bit errors (best scheme) = ' num2str(minBitErr)])
            
            figure(99), cla
            plot(1:trial, PrSignal(1:trial,iSNR), '-*'), hold on, grid on
            plot(1:trial, PrSignal_(1:trial,iSNR), '-o'), hold on, grid on
            plot(1:trial, ones(trial,1)*Np, 'k'), %overlay expected signal power
            xlabel('trial #'), ylabel('P_{signal}'), legend('Full CSI', 'Partial CSI')
            title(['Mean Rcvd Sig Pow: [' num2str(mean(PrSignal(1:trial)))...
                ', ' num2str(mean(PrSignal_(1:trial))) ']'])
            pause(0.1) %waiting a bit to get the plot window updated
        end
        
    end
    
    figure(18), cla
    semilogy(SNRvalues(1:iSNR), bitErrFullCSIR(1:iSNR), '-*'), hold on, grid on
    semilogy(SNRvalues(1:iSNR), bitErrPartCSIR(1:iSNR), '-o'), hold on, grid on
    xlabel('SNR (dB)'), ylabel('BER'),
    legend(['CP-OFDM (T = ' num2str(T) ', Full CSI)'],...
        ['CP-OFDM (T = ' num2str(T) ', Partial CSI)'])
    title(['\tau_{max} = ' num2str(tauMax)...
        ' ms, \alpha_{max} = '  num2str(alphaMax)...
        ', Np = ' num2str(Np)]);
end


%%
fc = 0;
H0 = delayScaleChanMatCPOFDM(hp, alphap,taup,N,L,fc,W,Fs);
