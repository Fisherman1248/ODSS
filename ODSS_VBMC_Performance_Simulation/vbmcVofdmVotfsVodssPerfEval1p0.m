% ------------------VBMC vs OFDM vs OTFS vs ODSS performance evaluation -------------------
%
% This code is written to assess the performance of the:
%
% (1) VBMC scheme that uses variable bandwidth chirps as subcarrier
% waveforms; the bandwidth of the waveform increases with subcarrier
% frequency.
%
% (2) OFDM scheme to use it as a benchmark for comparing the
% performance of the newly proposed VBMC scheme that uses variable
% bandwidth chirps as subcarrier waveforms; the bandwidth of the waveform
% increases with subcarrier frequency.
%
% (3) OTFS scheme as another benchmark for comparing the
% performance of the newly proposed VBMC scheme that uses variable
% bandwidth chirps as subcarrier waveforms; the bandwidth of the waveform
% increases with subcarrier frequency.
%
% A delay and scale spread channel with multiple paths is generated.
% Channel matrix, as seen at the receiver side after symbol extraction
% (through correlation rake-receiver), is computed using the analytical
% expression derived for VBMC, OFDM & OTFS Transmitter-to-Prop.
% Channel-to-Receiver system model where the propgation channel is a
% delay-scale channel
%
%Authors : (1) Arunkumar K. P.
%          (2) Chandra R. Murthy
%          (3) Muralikrishna P.
%Address : (1) Ph.D. Scholar,
%              Signal Processing for Communications Lab, ECE Department,
%              Indian Institute of Science, Bangalore, India-560 012.
%          (2) Professor,
%              Electrical Communication Engineering (ECE) Department,
%              Indian Institute of Science, Bangalore, India-560 012.
%          (3) Scientist,
%              Naval Physical Oceanographic Laboratory,
%              Kochi, India-682021.
% Email   : (1) arunkumar@iisc.ac.in
%
%
% Revision History
% Version : 1.1 is derived from vbmcVofdmVotfsPerfEval5px.m
% Created : 08-11-2022
%
% (1) ODSS modulation is now included for comparison
%
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% `````````````````````````````````````````````````````````````````````````
% Revision History of seed code (vbmcVofdmVotfsPerfEval5px.m)
% Version : 1.1.1 is derived from:
%           vbmcVofdmPerfParEval1p1 whose earliest verson was in turn
%           derived from
%           1.0 is derived from:
%           (1) vbmcPerfParEval1p2p1.m and
%           (2) ofdmPerfBenchParEval1p0p1.m
% Created : 21-12-2021
%
% (1) VBMC Code vbmcPerfParEval1p2p1 evolution history:
%            (v1.1) 09-01-2022: Modified chirp carriers
%            (v1.2) 13-01-2022: New channel matrix computing function -- to
%            speed up
%            14-01-2022: Found few issues, and their fix
%            15-01-2022: Modified the function delayScaleChanMatVBMC -- to
%            normalise the path modified subcarrier waveform dictionary Gp
%            s.t columns have unit norm => trace(Gp'*Gp) = number of
%            subcarriers
%            18-01-2022: Consistency in SNR definitions -- inband waveform
%            SNR is maintained when comparing across different modulation
%            schemes
%            (v2.0) 19-01-2022: parallel version using multiple cores
%            (v2.1) 19-01-2022: minor fixes -- averaging of some variables
%            left out in v2.0
% (2) OFDM code ofdmPerfBenchParEval1p0p1.m evolution history:
%            (v1.0)
%            18-01-2022: Consistency in SNR definitions -- inband waveform
%            SNR is maintained when comparing across different modulation
%            schemes
%            (v2.0)
%            19-01-2022: parallel version using multiple cores
%
% Revisions: (v1.1)
%            21-01-2022: saving fig and additional loop across of alphaMax
%            included
%            22-01-2022: loop across Np, tauMax included, ber and other
%            results are now saved, figures are also saved during the run
%            (v1.1.1)
%            22-01-2022: OTFS added
%            23-01-2022: OTFS partial CSIR based decoder made one-tap MMSE
%            equalizer; loop variable for Np, tauMax, alphaMax renamed to
%            fix the issue of replacement after one inner core looping
%            (v2.0)
%            26-01-2022: Created for modifying channel matrix computation,
%            Jain index evaluation for diagonal entries of the effective
%            channel matrix
%            29-01-2022: 
%            1) Simulations set up for symbol duration T = 80ms = betaB/dF.
%               For betaB = 2, dF = betaB/T = 1/40e-3 = 25 Hz
%            2) MMSE equalizer is now implemented at waveform level
%            (earlier version implemented MMSE equalizer at symbol level)
%            3) Reverted to MMSE Equalizer implemented at the symbol
%            measurement level, post MRC-MF -- since OTFS with single tap
%            (P-CSIR) equalizer was performing poorly (ber~0.5) when done
%            at waveform level! 
%            (v2.1.1.0) 
%            01-02-2022: Derived from vbmcVofdmVotfsPerfParEval2p1p1 for
%            evaluation on the personal macbook air at home. This version
%            is usually run after maskig parallel pool creation. Also, a
%            parallel pool size (poolsize) of at most 8 can be used on the
%            home macbook air. Usually poolsize is set to 2.
%            (v3.0)
%            01-02-2022: New version open to fix a major flaw! The length
%            of the VBMC waveform generated by the function vbmcWavesMod.m
%            exceeded the set duration by about 5 fold, compared to the
%            OFDM's transmitted waveform, which gave VBMC an SNR
%            advantage of 10*log10(5) ~ 7 dB at the MRC/MF output! The
%            mistake happened at line#18 of vbmcWavesMod.m which set: Tmax=
%            betaB/min(fHn-fLn) instead of Tmax = T. 
%            (v3.0) 
%            02-02-2022: parallel version for server runs derived from
%            vbmcVofdmVotfsPerfEval3p0.m
%            not so great results for VBMC
%            03-02-2022: to experiment the settings figured out on macbook
%            and setup extended runs on server
%            1) restored the waveform level MMSE equalization of symbols 
%            2) a MAJOR flaw found with the implementation of the one-tap
%            MMSE equalizer! MMSE equalizer matrix was constructed
%            improperly  without taking care of the fact that Matlab's
%            diag() function returns a vector (not matrix) whose entries
%            are the diagonal entries of the input matrix AND the  new
%            versions of Matlab allows the sum of a vector and matrix
%            whose number of rows match!!!
%
% INCORRECT usage: diag(A) will return a vector of diagonal entries of A,
% adding vector and matrix is a legitimate operation in newer versions of
% Matlab but will return a wrong result
%
% CORRECTED implementation: replace diag() with diag(diag())
%
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ VBMC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%                         One-tap MMSE equalizer
% INCORRECT implementation: yVpMMSEEq = TpV'*(diag(GrxV'*GrxV) + Pn0V*eye(N))\(GrxV'*rV(:)); %Equalizer MMSE One-Tap at the waveform level
% CORRECTED implementation: yVpMMSEEq = TpV'*(diag(diag(GrxV'*GrxV)) + Pn0V*eye(N))\(GrxV'*rV(:)); %Equalizer MMSE One-Tap at the waveform level                        
%
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ OFDM ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%                         One-tap MMSE equalizer
% INCORRECT implementation: yOpMMSEEq = TpO'*(diag(GrxO'*GrxO) + Pn0O*eye(Nfft))\(GrxO'*rO(:)); %One-Tap MMSE Equalizer at waveform level
% CORRECTED implementation: yOpMMSEEq = TpO'*(diag(diag(GrxO'*GrxO)) + Pn0O*eye(Nfft))\(GrxO'*rO(:)); %One-Tap MMSE Equalizer at waveform level
%
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ OTFS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%                         One-tap MMSE equalizer
% INCORRECT implementation: yTpMMSEEq = TpT'*(diag(GrxO'*GrxO) + Pn0T*eye(Nfft))\(GrxO'*rT(:)); %One-Tap MMSE Equalizer at waveform level
% CORRECTED implementation: yTpMMSEEq = (diag(diag((GrxO*TpT)'*(GrxO*TpT))) + Pn0T*eye(Nfft))\((GrxO*TpT)'*rT(:)); %One-Tap MMSE Equalizer at waveform level
%--------------------------------------------------------------------------
%            (vbmcVofdmVotfsICIandSIR.m)
%            04-02-2022: Included the computation of SIR. ICI, SIR and
%            Jain's index are now stored in the results.mat file
%            (vbmcVofdmVotfsPerfEval4x): "cross"-derived  from vbmcVofdmVotfsICIandSIR4x.m
%            06-02-2022: EESM based SIR-stat computation included for
%            several betaModn values, also min.,max., mean and median SIRs
%            are computed for comparing performance    
%            16-03-2022: Some betaM values od EESM computation altered: 5
%            and 30 included.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% This is a script under development.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% Setup the environment...
%  --plot paper settings, clear workspace, command window, parallel pool...
defpaper

clear
clc
close all

poolsize = 1;%56;
pOol = gcp('nocreate'); % If no pool, do not create new one.
if isempty(pOol)
    pOol = parpool(poolsize);
else    
    if pOol.NumWorkers ~= poolsize
        delete(gcp('nocreate'))
        pOol = parpool(poolsize);
    end
end
poolsize = pOol.NumWorkers;

%--------------------------------------------------------------------------
%% P.1) Channel parameters
NpList = 20;%20;%10;%[1 5 10 40];%20;%5;%10;%20; %number of paths from Tx to Rx
tauMaxList = 20e-3;%10e-3;%20e-3;%[20e-3 10e-3  40e-3 0 5e-3]; %max delay spread
alphaMaxList = 1.002;%[1.0 1.0005 1.002:0.001:1.005];%1.001;%[1.001 1.0 1.002:0.001:1.005 1.0005 1.0001];%:0.0001:1.0004]; %max scale spread

SNRlist = -9:3:21;%0:3:30; %in dB

%% P.2) System parameters
fL = 10e3;%start frequency, in Hz
fH = 20e3;%end frequency, in Hz
dF = 100;%100/2;%100;%bandwidth of the first subcarrier, in Hz
betaB = 2;%symbol extension factor
alphaMaxDes = 1.001; %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 1.001%max scale spread for which VBMC waveforms are designed

B = fH - fL; %bandwidth
FsV = 1.05*B;%1.040625*B;%1.0425*B;%sampling rate in VBMC system, in Hz
FsO = 1.05*B;%1.040625*B;%1.0425*B; %Sampling rate in OFDM system, in Hz

T = betaB/dF; %symbol duration, T > 1/dF
fc = (fH+fL)/2; %band center frequency

betaModn = [1.18 1.38 1.58 5.0 30.0];%[1.18 1.38 1.58 1.78 1.88];%beta parameters used to calculate EESM that depends on MCS -- trying out a few values

%% P.3) Monte-Carlo run parameters
Ntrials = 2*5000;%10000;%20000;%7000;%25000;% # of Monte Carlo trials to estimate BER
NBITERRMIN = 2*100;%200;%200; % # of bit errors to wait for, before closing the Monte Carlo run
NtrialMin = 2*100;%1000;%1000;%10000; %minimum number of trials to run (regardless of bit errors)
%--------------------------------------------------------------------------

%% Create VBMC subcarrier waveforms
% [ G, fLn, fHn ] = vbmcWaves( fL, fH, dF, T, FsV, alphaMaxDes );
% @@@@@@@@@@@@@@@@@@@@@@@@@@ [ G, fLn, fHn ] = vbmcWavesMod( fL, fH, dF, T, FsV, alphaMaxDes );
[ G, fLn, fHn ] = vbmcWavesModFix( fL, fH, dF, T, FsV, alphaMaxDes );
[M, N] = size(G);

SLL = 20;%77.1;%56.13;
gPuls = chebwin(M, SLL);
gPuls = [gPuls(2); gPuls(2:end-1); gPuls(end-1)];
% gPuls = ones(M,1);
GtxV = bsxfun(@times, G, gPuls);

gNorm = sqrt(trace(GtxV'*GtxV)/size(GtxV,2)); %normalizing factor to make trace(Gtx*Gtx') = size(Gtx,2)
GtxV = GtxV/gNorm;
Psubcar= diag(GtxV'*GtxV); %power of the subcarrier waveforms

ICIdB = ICI(abs(GtxV'*GtxV)); %transmit waveform ICI
SIRdB = SIR(abs(GtxV'*GtxV)); %transmit waveform SIR

figure(5), % Ambiguity function plot
imagesc(10*log10(abs(GtxV'*GtxV))), colorbar,
xlabel('subcarrier index (n)'), ylabel('subcarrier index (m)'),
title(['VBMC subcarrier correlation matrix. ICI_{max}='...
    num2str(round(10*ICIdB)/10) ' dB, SIR_{min}='...
    num2str(round(10*SIRdB)/10) ' dB.'] )

%% OFDM Waveform parameters
Nfft = floor(T*FsO);%betaB*N; %number of CP-OFDM subcarriers
% Lcp = ceil(FsO*tauMax); %CP-length
dataSubInds = 1 : Nfft/N : Nfft;%[1:betaB:N Nfft-N+1:betaB:Nfft];
if length(dataSubInds)~=N
    disp(['N = ' num2str(N)...
        'whereas length(dataSubInds) = ' num2str(length(dataSubInds))]);
    error('Number of subcarriers mounted with data symbols must match N (for fair comparison between VBMC/OFDM/OTFS)!')
end

%% Monte Carlo Trials
nTrialsPerBlock = poolsize;%~~~~~~~~~~~~~~~~~~~~2*poolsize;
disp(['Number of trials per trial-block : ' num2str(nTrialsPerBlock) ' ...'])

Ntrials = ceil(Ntrials/nTrialsPerBlock)*nTrialsPerBlock; %make Ntrials the least multiple of "nTrialsPerBlock" greater or equal to "Ntrials"
disp(...
    ['Number of trials (adjusted to the least multiple of nTrialsPerBlock >= original Ntrials) : '...
    num2str(Ntrials) ' ...']...
    );

nTrialBlocks = Ntrials/nTrialsPerBlock; %number of trial blocks to run
disp(['Number of trial-blocks : ' num2str(nTrialBlocks)])
% Ntrials = nTrialBlocks*nTrialsPerBlock

% VBMC performance monitors...
vbmcBERfullCSIRAvg = zeros(size(SNRlist)); %BER of the full complexity MMSE decoder using the entire effective channel matrix
vbmcBERpartCSIRAvg = zeros(size(SNRlist)); %BER of the one-tap MMSE decoder using only the diagonal entries of the effective channel matrix, D

vbmcPrSignalAvg = zeros(length(SNRlist),1);
vbmcPrSigCalAvg = zeros(length(SNRlist),1);
vbmcPySignalAvg = zeros(length(SNRlist),1);
vbmcPrNoiseAvg = zeros(length(SNRlist),1);
vbmcPyNoiseAvg = zeros(length(SNRlist),1);

vbmcMaxICIdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

vbmcSIRstatsdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList), 4+length(betaModn)+1);

vbmcSINRstatsdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList), 4+length(betaModn)+1);

vbmcSirJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

vbmcSinrJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

vbmcChanJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

% ODSS performance monitors...
odssBERfullCSIRAvg = zeros(size(SNRlist)); %BER of the full complexity MMSE decoder using the entire effective channel matrix
odssBERpartCSIRAvg = zeros(size(SNRlist)); %BER of the one-tap MMSE decoder using only the diagonal entries of the effective channel matrix, D

odssPrSignalAvg = zeros(length(SNRlist),1);
odssPrSigCalAvg = zeros(length(SNRlist),1);
odssPySignalAvg = zeros(length(SNRlist),1);
odssPrNoiseAvg = zeros(length(SNRlist),1);
odssPyNoiseAvg = zeros(length(SNRlist),1);

odssMaxICIdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

odssSIRstatsdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList), 4+length(betaModn)+1);

odssSINRstatsdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList), 4+length(betaModn)+1);

odssSirJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

odssSinrJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

odssChanJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

% OFDM performance monitors...
ofdmBERfullCSIRAvg = zeros(size(SNRlist)); %BER of the full complexity MMSE decoder using the entire effective channel matrix
ofdmBERpartCSIRAvg = zeros(size(SNRlist)); %BER of the one-tap MMSE decoder using only the diagonal entries of the effective channel matrix, D

ofdmPrSignalAvg = zeros(length(SNRlist),1);
ofdmPrSigCalAvg = zeros(length(SNRlist),1);
ofdmPySignalAvg = zeros(length(SNRlist),1);
ofdmPrNoiseAvg = zeros(length(SNRlist),1);
ofdmPyNoiseAvg = zeros(length(SNRlist),1);

ofdmMaxICIdBAvg = ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

ofdmSIRstatsdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList), 4+length(betaModn)+1);

ofdmSINRstatsdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList), 4+length(betaModn)+1);

ofdmSirJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

ofdmSinrJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

ofdmChanJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

% OTFS performance monitors...
otfsBERfullCSIRAvg = zeros(size(SNRlist)); %BER of the full complexity MMSE decoder using the entire effective channel matrix
otfsBERpartCSIRAvg = zeros(size(SNRlist)); %BER of the one-tap MMSE decoder using only the diagonal entries of the effective channel matrix, D

otfsPrSignalAvg = zeros(length(SNRlist),1);
otfsPrSigCalAvg = zeros(length(SNRlist),1);
otfsPySignalAvg = zeros(length(SNRlist),1);
otfsPrNoiseAvg = zeros(length(SNRlist),1);
otfsPyNoiseAvg = zeros(length(SNRlist),1);

otfsMaxICIdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

otfsSIRstatsdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList), 4+length(betaModn)+1);

otfsSINRstatsdBAvg =  ...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList), 4+length(betaModn)+1);

otfsSirJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

otfsSinrJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

otfsChanJainAvg =...
    zeros(length(NpList), length(tauMaxList), length(alphaMaxList));

% Setup figure panels
% -- to update plots in the background, use: set(0, 'CurrentFigure', figID)
%                                instead of: figure(figID)
berFig = 18; figure(berFig), pause(0.01); %BER plot
symbSnrFig = 27; figure(symbSnrFig), pause(0.01); %plot showing measured symbol SNR versus waveform SNR at the Rx for VBMC/OFDM
vbmcMeasFig = 99; figure(vbmcMeasFig), pause(0.01); %plots showing measured VBMC waveform, symbol and noise power, waveform/symbol SNR
ofdmMeasFig = 999; figure(ofdmMeasFig), pause(0.01); %plots showing measured OFDM waveform, symbol and noise power, waveform/symbol SNR
otfsMeasFig = 9999; figure(otfsMeasFig), pause(0.01); %plots showing measured OTFS waveform, symbol and noise power, waveform/symbol SNR
odssMeasFig = 99999; figure(odssMeasFig), pause(0.01); %plots showing measured ODSS waveform, symbol and noise power, waveform/symbol SNR
vbmcChanMatFig = 100; figure(vbmcChanMatFig), pause(0.01); % VBMC Rx Channel Matrix Image
ofdmChanMatFig = 1000; figure(ofdmChanMatFig), pause(0.01); % OFDM Rx Channel Matrix Image
otfsChanMatFig = 10000; figure(otfsChanMatFig), pause(0.01); % OTFS Rx Channel Matrix Image
odssChanMatFig = 100000; figure(odssChanMatFig), pause(0.01); % ODSS Rx Channel Matrix Image

displayRate = 30/60; %displays/plots are updated @ displayRate (minutes)
tDispl  = tic; % start the timer for displaying results
saveRate = 10; %results are saved @ saveRate (minutes)
tSave = tic; % start the timer for saving results
tExec = tic; % start timer for timing the loop execution

% Precoder matrix
TpV = eye(N);%exp(-1i*2*pi*(0:N-1)'*(0:N-1)/N)/sqrt(N);%eye(N);%
TpO = eye(Nfft);%exp(-1i*2*pi*(0:Nfft-1)'*(0:Nfft-1)/Nfft)/sqrt(Nfft);%eye(N);%eye(Nfft);%
TpT = exp(-1i*2*pi*(0:Nfft-1)'*(0:Nfft-1)/Nfft)/sqrt(Nfft);%eye(N);%eye(Nfft);%
TpD = exp(-1i*2*pi*(0:N-1)'*(0:N-1)/N)/sqrt(N);%eye(N);%eye(Nfft);%

% Runs...
iNp = 0;
for Np = NpList %loop across Np
    
    iNp = iNp + 1; iTau = 0; iAlpha = 0;
        
    disp(['Run for Np = ' num2str(Np) ' paths ...']);        
    
    for tauMax = tauMaxList %loop across tauMax
        
        iTau = iTau + 1; iAlpha = 0;
        
        disp(['Run for tauMax = ' num2str(tauMax*1e3) ' ms ...']);
        
        Lcp = ceil(FsO*tauMax); %CP-length                
        disp(['OFDM CP Length, Lcp = ' num2str(Lcp)]);
        
        % Inverse DFT matrix and its cyclic-prefixed version
        Fi = ifft(eye(Nfft));%inverse DFT matrix
        Fix = [Fi(end-Lcp+1:end,:);Fi]; %cycle prefixed inverse DFT matrix        
        
        for alphaMax = alphaMaxList %loop across alphaMax
            
            tExecAlpha = tic; %start timer for timing the execution of SNR loop for an alphaMax value
            
            iAlpha = iAlpha + 1;
            
            disp(['Run for alphaMax = ' num2str(alphaMax) ' ...']);
            
            for iSNR = 1:length(SNRlist) %loop across SNR
                
                tExecSnr = tic; %start timer for timing the Monte Carlo run execution for current SNR value
                
                SNRv = SNRlist(iSNR) - 10*log10(FsV/B); %read each inband SNR value, convert to SNR over the band
                SNRo = SNRlist(iSNR) - 10*log10(FsO/B); %read each inband SNR value, convert to SNR over the band
                
                vbmcPrSignal_ = zeros(nTrialBlocks, 1);
                vbmcPySignal_ = zeros(nTrialBlocks, 1);
                vbmcPrNoise_ = zeros(nTrialBlocks, 1);
                vbmcPyNoise_ = zeros(nTrialBlocks, 1);
                vbmcMaxICIdB_ = zeros(nTrialBlocks,1);
                vbmcSIRstatsdB_ = zeros(nTrialBlocks,4+length(betaModn)+1);
                vbmcSINRstatsdB_ = zeros(nTrialBlocks,4+length(betaModn)+1);
                vbmcSirJain_ = zeros(nTrialBlocks,1); 
                vbmcSinrJain_ = zeros(nTrialBlocks,1); 
                vbmcChanJain_ = zeros(nTrialBlocks,1);
                
                odssPrSignal_ = zeros(nTrialBlocks, 1);
                odssPySignal_ = zeros(nTrialBlocks, 1);
                odssPrNoise_ = zeros(nTrialBlocks, 1);
                odssPyNoise_ = zeros(nTrialBlocks, 1);
                odssMaxICIdB_ = zeros(nTrialBlocks,1);
                odssSIRstatsdB_ = zeros(nTrialBlocks,4+length(betaModn)+1);
                odssSINRstatsdB_ = zeros(nTrialBlocks,4+length(betaModn)+1);
                odssSirJain_ = zeros(nTrialBlocks,1); 
                odssSinrJain_ = zeros(nTrialBlocks,1); 
                odssChanJain_ = zeros(nTrialBlocks,1);                
                
                ofdmPrSignal_ = zeros(nTrialBlocks, 1);
                ofdmPySignal_ = zeros(nTrialBlocks, 1);
                ofdmPrNoise_ = zeros(nTrialBlocks, 1);
                ofdmPyNoise_ = zeros(nTrialBlocks, 1);
                ofdmMaxICIdB_ = zeros(nTrialBlocks,1);
                ofdmSIRstatsdB_ = zeros(nTrialBlocks,4+length(betaModn)+1);
                ofdmSINRstatsdB_ = zeros(nTrialBlocks,4+length(betaModn)+1);
                ofdmSirJain_ = zeros(nTrialBlocks,1); 
                ofdmSinrJain_ = zeros(nTrialBlocks,1);                 
                ofdmChanJain_ = zeros(nTrialBlocks,1); 
                
                otfsPrSignal_ = zeros(nTrialBlocks, 1);
                otfsPySignal_ = zeros(nTrialBlocks, 1);
                otfsPrNoise_ = zeros(nTrialBlocks, 1);
                otfsPyNoise_ = zeros(nTrialBlocks, 1);
                otfsMaxICIdB_ = zeros(nTrialBlocks,1);    
                otfsSIRstatsdB_ = zeros(nTrialBlocks,4+length(betaModn)+1);
                otfsSINRstatsdB_ = zeros(nTrialBlocks,4+length(betaModn)+1);
                otfsSirJain_ = zeros(nTrialBlocks,1); 
                otfsSinrJain_ = zeros(nTrialBlocks,1); 
                otfsChanJain_ = zeros(nTrialBlocks,1); 
                
                for iTrialBlock = 1 : nTrialBlocks %run each trial-block
                    
                    tExecTrialBlock = tic; %start timer for timing the Monte Carlo run execution for current SNR value
                    
                    %------------------------------------------------------------------
                    vbmcBERfullCSIR = zeros(nTrialsPerBlock, 1); %BER in each trial run of the full complexity MMSE decoder
                    vbmcBERpartCSIR = zeros(nTrialsPerBlock, 1); %BER in each trial run of the one-tap MMSE decoder
                    
                    vbmcPrSignal = zeros(nTrialsPerBlock, 1);
                    vbmcPrSigCal = zeros(nTrialsPerBlock, 1);
                    vbmcPySignal = zeros(nTrialsPerBlock, 1);
                    vbmcPrNoise = zeros(nTrialsPerBlock, 1);
                    vbmcPyNoise = zeros(nTrialsPerBlock, 1);
                    
                    vbmcMaxICIdB = zeros(nTrialsPerBlock,1);
                    vbmcSIRstatsdB = zeros(nTrialsPerBlock,4+length(betaModn)+1);
                    vbmcSINRstatsdB = zeros(nTrialsPerBlock,4+length(betaModn)+1);                    
                    vbmcSirJain = zeros(nTrialsPerBlock,1);
                    vbmcSinrJain = zeros(nTrialsPerBlock,1);                    
                    vbmcChanJain = zeros(nTrialsPerBlock,1);
                    %------------------------------------------------------------------
                    odssBERfullCSIR = zeros(nTrialsPerBlock, 1); %BER in each trial run of the full complexity MMSE decoder
                    odssBERpartCSIR = zeros(nTrialsPerBlock, 1); %BER in each trial run of the one-tap MMSE decoder
                    
                    odssPrSignal = zeros(nTrialsPerBlock, 1);
                    odssPrSigCal = zeros(nTrialsPerBlock, 1);
                    odssPySignal = zeros(nTrialsPerBlock, 1);
                    odssPrNoise = zeros(nTrialsPerBlock, 1);
                    odssPyNoise = zeros(nTrialsPerBlock, 1);
                    
                    odssMaxICIdB = zeros(nTrialsPerBlock,1);
                    odssSIRstatsdB = zeros(nTrialsPerBlock,4+length(betaModn)+1);
                    odssSINRstatsdB = zeros(nTrialsPerBlock,4+length(betaModn)+1);                    
                    odssSirJain = zeros(nTrialsPerBlock,1);
                    odssSinrJain = zeros(nTrialsPerBlock,1);                    
                    odssChanJain = zeros(nTrialsPerBlock,1);                    
                    %------------------------------------------------------------------
                    ofdmBERfullCSIR = zeros(nTrialsPerBlock, 1); %BER in each trial run of the full complexity MMSE decoder
                    ofdmBERpartCSIR = zeros(nTrialsPerBlock, 1); %BER in each trial run of the one-tap MMSE decoder
                    
                    ofdmPrSignal = zeros(nTrialsPerBlock, 1);
                    ofdmPrSigCal = zeros(nTrialsPerBlock, 1);
                    ofdmPySignal = zeros(nTrialsPerBlock, 1);
                    ofdmPrNoise = zeros(nTrialsPerBlock, 1);
                    ofdmPyNoise = zeros(nTrialsPerBlock, 1);
                    
                    ofdmMaxICIdB = zeros(nTrialsPerBlock,1);
                    ofdmSIRstatsdB = zeros(nTrialsPerBlock,4+length(betaModn)+1);
                    ofdmSINRstatsdB = zeros(nTrialsPerBlock,4+length(betaModn)+1);
                    ofdmSirJain = zeros(nTrialsPerBlock,1);
                    ofdmSinrJain = zeros(nTrialsPerBlock,1);                    
                    ofdmChanJain = zeros(nTrialsPerBlock,1);
                    %------------------------------------------------------------------
                    otfsBERfullCSIR = zeros(nTrialsPerBlock, 1); %BER in each trial run of the full complexity MMSE decoder
                    otfsBERpartCSIR = zeros(nTrialsPerBlock, 1); %BER in each trial run of the one-tap MMSE decoder
                    
                    otfsPrSignal = zeros(nTrialsPerBlock, 1);
                    otfsPrSigCal = zeros(nTrialsPerBlock, 1);
                    otfsPySignal = zeros(nTrialsPerBlock, 1);
                    otfsPrNoise = zeros(nTrialsPerBlock, 1);
                    otfsPyNoise = zeros(nTrialsPerBlock, 1);
                    
                    otfsMaxICIdB = zeros(nTrialsPerBlock,1);
                    otfsSIRstatsdB = zeros(nTrialsPerBlock,4+length(betaModn)+1);
                    otfsSINRstatsdB = zeros(nTrialsPerBlock,4+length(betaModn)+1);                    
                    otfsSirJain = zeros(nTrialsPerBlock,1);
                    otfsSinrJain = zeros(nTrialsPerBlock,1);
                    otfsChanJain = zeros(nTrialsPerBlock,1);
                    %------------------------------------------------------------------
                    PnoisetV = zeros(nTrialsPerBlock,1);
                    HtV = zeros(N,N,nTrialsPerBlock);
                    Ps0V = Np*N; %true or expected signal (waveform) power at receiver
                    %------------------------------------------------------------------
                    PnoisetD = zeros(nTrialsPerBlock,1);
                    HtD = zeros(N,N,nTrialsPerBlock);
                    Ps0D = Np*N; %true or expected signal (waveform) power at receiver                    
                    %------------------------------------------------------------------
                    PnoisetO = zeros(nTrialsPerBlock,1);
                    HtO = zeros(N,N,nTrialsPerBlock);
                    Ps0O = Np * (40*sqrt(3)*N/70);%true or expected signal (waveform) power at receiver
                    %------------------------------------------------------------------
                    PnoisetT = zeros(nTrialsPerBlock,1);
                    HtT = zeros(N,N,nTrialsPerBlock);
                    Ps0T = Np * (40*sqrt(3)*N/70);%true or expected signal (waveform) power at receiver                                        
                    %------------------------------------------------------------------
                    
                    disp(['SNR = ' num2str(SNRlist(iSNR)) ' dB,'...
                        ' alphaMax = ' num2str(alphaMax) ','...
                        ' tauMax = ' num2str(tauMax*1e3) ' ms,'...
                        ' Np = ' num2str(Np) ' paths,'...
                        ' trial-block = ' num2str(iTrialBlock)...
                        ' of ' num2str(nTrialBlocks) ' trial-blocks.']);
                    
                    parfor trial = 1:nTrialsPerBlock                    
%                         disp(['Trial ' num2str(trial) '/'...
%                             num2str(nTrialsPerBlock) ' in block '...
%                             num2str(iTrialBlock) '/' num2str(nTrialBlocks) ':'...
%                             ' SNR = ' num2str(SNRlist(iSNR)) ' dB,'...
%                             ' alphaMax = ' num2str(alphaMax) ','...
%                             ' tauMax = ' num2str(tauMax*1e3) ' ms,'...
%                             ' Np = ' num2str(Np) ' paths'...
%                             ]);
                        
                        %% Transmit Side: Modulator (mounts data symbols)
                        x = (randi(2,[N 1])-1-1/2)*2; %symbols/bits to transmit through VBMC/OFDM
                        
                        %^^^^^^^^^^^^^^^^^^^^ VBMC ^^^^^^^^^^^^^^^^^^^^^^^^
                        xV = x; %transfer bits to be transmitted to VBMC symbol vector
                        s = GtxV*(TpV*xV);%transmitted signal -- VBMC modulator output
                        %^^^^^^^^^^^^^^^^^^^^ ODSS ^^^^^^^^^^^^^^^^^^^^^^^^
                        xD = x; %transfer bits to be transmitted to ODSS symbol vector
                        sD = GtxV*(TpD*xD);%transmitted signal -- ODSS modulator output                        
                        %^^^^^^^^^^^^^^^^^^^^ OFDM ^^^^^^^^^^^^^^^^^^^^^^^^                        
                        xO = zeros(Nfft,1);
                        xO(dataSubInds) = x;%transfer bits to be transmitted to OFDM symbol vector
                        %^^^^^^^^^^^^^^^^^^^^ OTFS ^^^^^^^^^^^^^^^^^^^^^^^^                        
                        xT = zeros(Nfft,1);
                        xT(dataSubInds) = x;%transfer bits to be transmitted to OFDM symbol vector
                        %--------------------------------------------------                        
                        
                        %         figure(6), clf,
                        %         subplot(211), stem(x), grid on,
                        %         xlabel('symbol index'), ylabel('symbol value'), title('transmitted symbols')
                        %         subplot(212), plot(t, real(s)), grid on, hold on, plot(t, imag(s))
                        %         xlabel('t (s)'), ylabel('s(t)'), title('transmitted signal')
                        %--------------------------------------------------
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
                        %--------------------------------------------------
                        %% Rx waveform level channel matrix
                        [Hsinc] = delayScaleChannelBasis(hp, taup, alphap,...
                            FsO, fc, FsO, (Nfft+Lcp), tauMax); % delay-scale channel matrix, H
                        
                        %^^^^^^^^^^^^^^^^^^^^ VBMC ^^^^^^^^^^^^^^^^^^^^^^^^
%                         GrxV = delayScaleChanMatVBMCmodFix(hp, alphap, taup, T, fLn, fHn, betaB, FsV, tauMax);%, SLL);
%                         GrxV = delayScaleChanMatVBMC(hp, alphap, taup, T, fLn, fHn, FsV, tauMax);%, SLL);
                        GrxV = Hsinc*[GtxV;zeros(Lcp,N)];%Grx = H*Gtx

                        %^^^^^^^^^^^^^^^^^^^^ ODSS ^^^^^^^^^^^^^^^^^^^^^^^^
%                         GrxV = delayScaleChanMatVBMCmodFix(hp, alphap, taup, T, fLn, fHn, betaB, FsV, tauMax);%, SLL);
%                         GrxV = delayScaleChanMatVBMC(hp, alphap, taup, T, fLn, fHn, FsV, tauMax);%, SLL);
                        GrxD = Hsinc*[GtxV;zeros(Lcp,N)];%Grx = H*Gtx
                        
                        %^^^^^^^^^^^^^^^^^^^ OFDM/OTFS ^^^^^^^^^^^^^^^^^^^^
%                         GrxO = delayScaleChanMatCPOFDM(hp, alphap, taup, Nfft, Lcp, fc, B, FsO);                        
                        GrxO = fft(Hsinc(Lcp+1:Lcp+Nfft, 1:Lcp+Nfft)*Fix); %select only the samples corresponding to Nfft symbols
                                                                        
                        %% Received Signal
                        %--------------------------------------------------------------
                        rsV = GrxV*(TpV*xV);%waveform reaching the receiver after propagation from the transmitter 
                        %             Ps0V = Np*N; %true or expected signal (waveform) power at receiver
                        
                        Pn0V = Ps0V*10^(-SNRv/10)/length(rsV); %noise power (at waveform level) for the specified SNR
                        wV = sqrt(Pn0V)*(randn(size(rsV))+1i*randn(size(rsV)))/sqrt(2); %AWGN
                        
                        rV = rsV + wV; %received signal after (additive) noise corruption
                        
                        PsV = rsV'*rsV;%estimated received signal power (at waveform level) 
                        vbmcPrSignal(trial) = PsV; %record the measured signal power
                        
                        PnV = wV'*wV; %estimated receiver noise power (at waveform level) 
                        vbmcPrNoise(trial) = PnV;%record the measured noise power
                        
                        vbmcPrSigCal(trial) = trace(GrxV'*GrxV);%alternate evaluation

                        %--------------------------------------------------------------
                        rsD = GrxD*(TpD*xD);%waveform reaching the receiver after propagation from the transmitter 
                        %             Ps0D = Np*N; %true or expected signal (waveform) power at receiver
                        
                        Pn0D = Ps0V*10^(-SNRv/10)/length(rsD); %noise power (at waveform level) for the specified SNR
                        wD = sqrt(Pn0D)*(randn(size(rsD))+1i*randn(size(rsD)))/sqrt(2); %AWGN
                        
                        rD = rsD + wD; %received signal after (additive) noise corruption
                        
                        PsD = rsD'*rsD;%estimated received signal power (at waveform level) 
                        odssPrSignal(trial) = PsD; %record the measured signal power
                        
                        PnD = wD'*wD; %estimated receiver noise power (at waveform level) 
                        odssPrNoise(trial) = PnD;%record the measured noise power
                        
                        odssPrSigCal(trial) = trace(GrxD'*GrxD);%alternate evaluation
                        
                        %--------------------------------------------------------------
                        rsO = GrxO*(TpO*xO); %received waveform
                        %         Ps0 = Np * (40*sqrt(3)*N/70);%true or expected signal (waveform) power at receiver
                        
                        Pn0O = Ps0O*10^(-SNRo/10)/Nfft;%noise (waveform-level) power for the specified SNR
                        wO = sqrt(Pn0O)*(randn(size(rsO))+1i*randn(size(rsO)))/sqrt(2); %AWGN
                        
                        rO = rsO + wO; %received waveform after (additive) noise corruption
                        
                        PsO = rsO'*rsO; %estimated (empirical) received signal (waveform level) power
                        ofdmPrSignal(trial) = PsO;  %monitor/log received signal power measurement
                        
                        PnO = wO'*wO;   %estimated (empirical) receiver noise power
                        ofdmPrNoise(trial) = PnO;  %monitor/log received signal power measurement
                        
                        ofdmPrSigCal(trial) = trace(GrxO'*GrxO)/betaB;%alternate evaluation
                        %--------------------------------------------------------------
                        %--------------------------------------------------------------
                        rsT = GrxO*(TpT*xT); %received waveform
                        %         Ps0 = Np * (40*sqrt(3)*N/70);%true or expected signal (waveform) power at receiver
                        
                        Pn0T = Ps0T*10^(-SNRo/10)/Nfft;%noise (waveform-level) power for the specified SNR
                        wT = sqrt(Pn0T)*(randn(size(rsT))+1i*randn(size(rsT)))/sqrt(2); %AWGN
                        
                        rT = rsT + wT; %received waveform after (additive) noise corruption
                        
                        PsT = rsT'*rsT; %estimated (empirical) received signal (waveform level) power
                        otfsPrSignal(trial) = PsT;  %monitor/log received signal power measurement
                        
                        PnT = wT'*wT;   %estimated (empirical) receiver noise power
                        otfsPrNoise(trial) = PnT;  %monitor/log received signal power measurement
                        
                        otfsPrSigCal(trial) = trace(GrxO'*GrxO)/betaB;%alternate evaluation
                        %--------------------------------------------------------------
                        
                        %         figure(9), clf,%plots the transmitted signal and the received signal
                        %         tr = (0:length(r)-1)/Fs; %sampling times for the received signal
                        %         plot(tr, real(r)), grid on, hold on, plot(tr, imag(r))
                        %         xlabel('t (s)'), ylabel('r(t)'), legend('Re', 'Im')
                        %         title('Received signal after passing s(t) through delay-scale channel')
                        
                        %% Receiver processing
                        %----------------------------------------------------------
                        yV = GrxV'*rV(:); %VBMC symbol measurements, post MRC+MF
                        
                        vbmcPyNoise(trial) = (GrxV'*wV)'*(GrxV'*wV); %noise power after matched filtering
                        vbmcPySignal(trial) = (GrxV'*rsV(:))'*(GrxV'*rsV(:)); %signal power after matched filtering
                        
                        PnoiseV = Ps0V*Pn0V/N;%norm(y(:)-H*(Tp*x(:)),2)^2/N;% Pnoise = E{trace(Grx'*Grx)}*Pn0/N = Ps0*Pn0/N;% mean(PyNoise(1:trial,iSNR))/N
                        PnoisetV(trial) = PnoiseV;                        
                        
                        Hv = GrxV'*GrxV; %effective channel matrix, as seen after symbol measurements
                        vbmcMaxICIdB(trial) = ICI(sqrt(abs(Hv)));%maxnetICI(sqrt(abs(Hv)));%~~~~~~~ICI(sqrt(abs(Hv)));
                        [vbmcSIRstatsdB(trial,:), vbmcSirJain(trial)] =...
                            SIRstats(sqrt(abs(Hv)),betaModn);
                        [vbmcSINRstatsdB(trial,:), vbmcSinrJain(trial)] =...
                            SINRstats(sqrt(abs(Hv)),PnoiseV,betaModn);
                        
                        HtV(:,:,trial) = Hv./diag(Hv);
                        
                        %----------------------------------------------------------
                        yD = GrxD'*rD(:); %ODSS symbol measurements, post MRC+MF
                        
                        odssPyNoise(trial) = (GrxD'*wD)'*(GrxD'*wD); %noise power after matched filtering
                        odssPySignal(trial) = (GrxD'*rsD(:))'*(GrxD'*rsD(:)); %signal power after matched filtering
                        
                        PnoiseD = Ps0D*Pn0D/N;%norm(y(:)-H*(Tp*x(:)),2)^2/N;% Pnoise = E{trace(Grx'*Grx)}*Pn0/N = Ps0*Pn0/N;% mean(PyNoise(1:trial,iSNR))/N
                        PnoisetD(trial) = PnoiseD;                        
                        
                        Hd = GrxD'*GrxD; %effective channel matrix, as seen after symbol measurements
                        odssMaxICIdB(trial) = ICI(sqrt(abs(Hd)));%maxnetICI(sqrt(abs(Hv)));%~~~~~~~ICI(sqrt(abs(Hv)));
                        [odssSIRstatsdB(trial,:), odssSirJain(trial)] =...
                            SIRstats(sqrt(abs(Hd)),betaModn);
                        [odssSINRstatsdB(trial,:), odssSinrJain(trial)] =...
                            SINRstats(sqrt(abs(Hd)),PnoiseD,betaModn);
                        
                        HtD(:,:,trial) = Hd./diag(Hd);                        
                        
                        %----------------------------------------------------------
                        yO = GrxO'*rO(:); %OFDM symbol measurements, post MRC+MF
                        
                        ofdmPyNoise(trial) = (GrxO'*wO)'*(GrxO'*wO); %noise power after matched filtering
                        ofdmPySignal(trial) = (GrxO'*rsO(:))'*(GrxO'*rsO(:)); %signal power after matched filtering
                        
                        PnoiseO = Ps0O*betaB*Pn0O/N*sqrt(2);%norm(y(:)-H*(Tp*x(:)),2)^2/N;% Pnoise = E{trace(Grx'*Grx)}*Pn0/N = Ps0*betaB*Pn0/N; % mean(PyNoise(1:trial,iSNR))/N
                        PnoisetO(trial) = PnoiseO;                        
                        
                        HO = GrxO'*GrxO;%effective channel matrix, as seen after symbol measurements
                        ofdmMaxICIdB(trial) = ICI(sqrt(abs(HO)));%maxICI(sqrt(abs(HO)));%~~~~~~~~~~~ICI(sqrt(abs(HO)));
                        [ofdmSIRstatsdB(trial,:), ofdmSirJain(trial)] =...
                            SIRstats(sqrt(abs(HO)),betaModn);
                        [ofdmSINRstatsdB(trial,:), ofdmSinrJain(trial)] =...
                            SINRstats(sqrt(abs(HO)),PnoiseO,betaModn);
                        
                        HtO(:,:,trial) = HO(dataSubInds,dataSubInds)./...
                            diag(HO(dataSubInds,dataSubInds));
                        
                        %----------------------------------------------------------
                        yT = (GrxO*TpT)'*rT(:); %OTFS symbol measurements, post MRC+MF
                        
                        otfsPyNoise(trial) = ((GrxO*TpT)'*wT)'*((GrxO*TpT)'*wT); %noise power after matched filtering
                        otfsPySignal(trial) = ((GrxO*TpT)'*rsT(:))'*((GrxO*TpT)'*rsT(:)); %signal power after matched filtering
                        
                        PnoiseT = Ps0T*betaB*Pn0T/N*sqrt(2);%norm(y(:)-H*(Tp*x(:)),2)^2/N;% Pnoise = E{trace(Grx'*Grx)}*Pn0/N = Ps0*betaB*Pn0/N; % mean(PyNoise(1:trial,iSNR))/N
                        PnoisetT(trial) = PnoiseT;                        
                        
                        HT = (GrxO*TpT)'*(GrxO*TpT);%effective channel matrix, as seen after symbol measurements
                        otfsMaxICIdB(trial) = ICI(sqrt(abs(HT)));%maxnetICI(sqrt(abs(HO)));%~~~~~~~~~~~~~~~~~ICI(sqrt(abs(HO)));
                        [otfsSIRstatsdB(trial,:), otfsSirJain(trial)] =...
                            SIRstats(sqrt(abs(HT)),betaModn);
                        [otfsSINRstatsdB(trial,:), otfsSinrJain(trial)] =...
                            SINRstats(sqrt(abs(HT)),PnoiseT,betaModn);
                        
                        HtT(:,:,trial) = HT(dataSubInds,dataSubInds)./...
                            diag(HT(dataSubInds,dataSubInds));
                        
                        %--------------------------------------------------
                        
                        %% Jain Index of the Effective Channel Matrix of each modulation scheme
                        vbmcChanJain(trial) = JainIndex(diag(abs(Hv*TpV)));
                        odssChanJain(trial) = JainIndex(diag(abs(Hd*TpD)));
                        ofdmChanJain(trial) = JainIndex(diag(abs(HO*TpO)));
                        otfsChanJain(trial) = JainIndex(diag(abs(HO*TpT)));                        
                        
                        %% (1) Symbol recovery using full CSIR
                        %--------------------------------------------------
                        %         yfZFEq = (H*Tp)\y;%(diag(diag(Hmat))*(Gtx*Gtx'))\Dr;%ZF Equalizer
                        yVfMMSEEq = ((GrxV*TpV)'*(GrxV*TpV) + Pn0V*eye(N))\((GrxV*TpV)'*rV(:)); %MMSE Equalizer at the waveform level
%                         yVfMMSEEq = ((Hv*TpV)'*(Hv*TpV) + PnoiseV*eye(N))\((Hv*TpV)'*yV(:)); %MMSE Equalizer at the symbol level                        
                        
                        vbmcBERfullCSIR(trial) = sum(ne(xV>0, real(yVfMMSEEq)>0))/length(xV);
                        %--------------------------------------------------
                        %         yfZFEq = (H*Tp)\y;%(diag(diag(Hmat))*(Gtx*Gtx'))\Dr;%ZF Equalizer
                        yDfMMSEEq = ((GrxD*TpD)'*(GrxD*TpD) + Pn0D*eye(N))\((GrxD*TpD)'*rD(:)); %MMSE Equalizer at the waveform level
%                         yDfMMSEEq = ((Hd*TpD)'*(Hd*TpD) + PnoiseD*eye(N))\((Hd*TpD)'*yD(:)); %MMSE Equalizer at the symbol level                        
                        
                        odssBERfullCSIR(trial) = sum(ne(xD>0, real(yDfMMSEEq)>0))/length(xD);                        
                        %--------------------------------------------------
                        %         yfZFEq = (H*Tp)\y;%(diag(diag(Hmat))*(Gtx*Gtx'))\Dr;%LS Decoder
                        yOfMMSEEq = ((GrxO*TpO)'*(GrxO*TpO) + Pn0O*eye(Nfft))\((GrxO*TpO)'*rO(:)); %MMSE Equalizer at the waveform level
%                         yOfMMSEEq = ((HO*TpO)'*(HO*TpO) + PnoiseO*eye(Nfft))\((HO*TpO)'*yO(:)); %MMSE Equalizer at the symbol level
                        
                        ofdmBERfullCSIR(trial) = sum(ne(xO(dataSubInds)>0,...
                            real(yOfMMSEEq(dataSubInds))>0)) / length(xO(dataSubInds));
                        %--------------------------------------------------
                        %         yfZFEq = (H*Tp)\y;%(diag(diag(Hmat))*(Gtx*Gtx'))\Dr;%LS Decoder
                        yTfMMSEEq = ((GrxO*TpT)'*(GrxO*TpT) + Pn0T*eye(Nfft))\((GrxO*TpT)'*rT(:)); %MMSE Equalizer at the waveform level
%                         yTfMMSEEq = ((HO*TpT)'*(HO*TpT) + PnoiseT*eye(Nfft))\((HO*TpT)'*yT(:)); %MMSE Equalizer at the symbol level
                        
                        otfsBERfullCSIR(trial) = sum(ne(xT(dataSubInds)>0,...
                            real(yTfMMSEEq(dataSubInds))>0)) / length(xT(dataSubInds));
                        %----------------------------------------------------------
                        
                        %         figure(12),
                        %         subplot(211), cla, stem(x), hold on, stem(real(DrfEq),'*'), grid on, stem(imag(DrfEq),'.') % stem(abs(DrfEq_),'*'), grid on
                        %         xlabel('bit index'), ylabel('bit value'), ylim([-2 2])
                        %         legend('Transmitted bit', 'Received (soft) bit: Real', 'Received (soft) bit: Imag')
                        %         title('Symbol recovery using full CSI')   
                        
                        %--------------------------------------------------                        
                        %% (3) Symbol recovery using partial CSIR (only diagonal entry of Hmat)     
                        
                        % ^^^^^^^^^^^^^^^^^^^^^^VBMC^^^^^^^^^^^^^^^^^^^^^^^                        
                        %         ypZFEq = Tp\(y./diag(H)); %ZF Equalizer
%@@@@@@@@@@@@@@@@@@@@@@ yVpMMSEEq = TpV'*(diag(GrxV'*GrxV) + Pn0V*eye(N))\(GrxV'*rV(:)); %Equalizer MMSE One-Tap at the waveform level
                        yVpMMSEEq = (diag(diag((GrxV*TpV)'*(GrxV*TpV))) + Pn0V*eye(N))\((GrxV*TpV)'*rV(:)); %Equalizer MMSE One-Tap at the waveform level
                        
%                         yVpMMSEEq = ((diag(diag(Hv))*TpV)'*(diag(diag(Hv))*TpV) +...
%                             PnoiseV*eye(N))\((diag(diag(Hv))*TpV)'*yV(:)); %One-Tap MMSE Equalizer at symbol level                        
                        
                        vbmcBERpartCSIR(trial) = sum(ne(xV>0, real(yVpMMSEEq)>0))/length(xV);

                        % ^^^^^^^^^^^^^^^^^^^^^^ODSS^^^^^^^^^^^^^^^^^^^^^^^                        
                        %         ypZFEq = Tp\(y./diag(H)); %ZF Equalizer
%@@@@@@@@@@@@@@@@@@@@@@ yDpMMSEEq = TpD'*(diag(GrxD'*GrxD) + Pn0D*eye(N))\(Grxd'*rD(:)); %Equalizer MMSE One-Tap at the waveform level
%                         yDpMMSEEq = (diag(diag((GrxD*TpD)'*(GrxD*TpD))) + Pn0D*eye(N))\((GrxD*TpD)'*rD(:)); %Equalizer MMSE One-Tap at the waveform level
                        yDpMMSEEq = TpD'*((diag(diag((GrxD)'*(GrxD))) + Pn0D*eye(N))\((GrxD)'*rD(:))); %Equalizer MMSE One-Tap at the waveform level
                        
%                         yVpMMSEEq = ((diag(diag(Hv))*TpV)'*(diag(diag(Hv))*TpV) +...
%                             PnoiseV*eye(N))\((diag(diag(Hv))*TpV)'*yV(:)); %One-Tap MMSE Equalizer at symbol level                        
                        
                        odssBERpartCSIR(trial) = sum(ne(xD>0, real(yDpMMSEEq)>0))/length(xD);
                        
                        % ^^^^^^^^^^^^^^^^^^^^^ OFDM ^^^^^^^^^^^^^^^^^^^^^^
                        %         ypZFEq = Tp\(y./diag(H)); %ZF Equalizer
%@@@@@@@@@@@@@@@@@@@@@@ yOpMMSEEq = TpO'*(diag(GrxO'*GrxO) + Pn0O*eye(Nfft))\(GrxO'*rO(:)); %One-Tap MMSE Equalizer at waveform level
                        yOpMMSEEq = (diag(diag((GrxO*TpO)'*(GrxO*TpO))) + Pn0O*eye(Nfft))\((GrxO*TpO)'*rO(:)); %One-Tap MMSE Equalizer at waveform level
                        
%                         yOpMMSEEq = ((diag(diag(HO))*TpO)'*(diag(diag(HO))*TpO) +...
%                             PnoiseO*eye(Nfft))\((diag(diag(HO))*TpO)'*yO(:)); %One-Tap MMSE Equalizer at symbol level                        
                        
                        ofdmBERpartCSIR(trial) =  sum(ne(xO(dataSubInds)>0,...
                            real(yOpMMSEEq(dataSubInds))>0)) / length(xO(dataSubInds));
                        
                        % ^^^^^^^^^^^^^^^^^^^^^ OTFS ^^^^^^^^^^^^^^^^^^^^^^
                        %         ypZFEq = Tp\(y./diag(H)); %ZF Equalizer
%@@@@@@@@@@@@@@@@@@@@@@ yTpMMSEEq = TpT'*(diag(GrxO'*GrxO) + Pn0T*eye(Nfft))\(GrxO'*rT(:)); %One-Tap MMSE Equalizer at waveform level
%                         yTpMMSEEq = (diag(diag((GrxO*TpT)'*(GrxO*TpT))) + Pn0T*eye(Nfft))\((GrxO*TpT)'*rT(:)); %One-Tap MMSE Equalizer at waveform level
                        yTpMMSEEq = TpT'*((diag(diag((GrxO)'*(GrxO))) + Pn0T*eye(Nfft))\((GrxO)'*rT(:))); %One-Tap MMSE Equalizer at waveform level                        
                        
%                        %##################################################
%                         yTpMMSEEq = ((diag(diag(HO))*TpT)'*(diag(diag(HO))*TpT) +...
%                             PnoiseT*eye(Nfft))\((diag(diag(HO))*TpT)'*yT(:));
%                             %MMSE Equalizer using diagonal approximation of
%                             %the effective channel matrix HO - this
%                             %equalizer is NOT ONE-TAP though!
%                        %##################################################

%                         yTpMMSEEq = TpT'*( ((diag(diag(HO))'*diag(diag(HO))) +...
%                             PnoiseT*eye(Nfft))\((diag(diag(HO)))'*yT(:))); %One-Tap MMSE Equalizer at symbol level

                        otfsBERpartCSIR(trial) =  sum(ne(xT(dataSubInds)>0,...
                            real(yTpMMSEEq(dataSubInds))>0)) / length(xT(dataSubInds));
                        %--------------------------------------------------
                        
                        %         figure(12),
                        %         subplot(212), cla, stem(x), hold on, stem(real(DrpEq),'*'), grid on, stem(imag(DrpEq),'.') %stem(abs(DrpEq),'*'), grid on
                        %         xlabel('bit index'), ylabel('bit value'), ylim([-2 2])
                        %         legend('Transmitted bit', 'Received (soft) bit: Real', 'Received (soft) bit: Imag')
                        %         title('Symbol recovery using single tap equalizer')
                        %         pause(0.01);
                        
                        %         if (bitErrPartCSIR(iSNR)>0) && (SNR == inf)
                        %             disp('Bit Error found despite SNR = infty! Terminating execution to leave way for analysing...')
                        %             break;
                        %         end
                    end
                    
                    toc(tExecTrialBlock),
                    
                    %% Compute all averages till the most recent trial
                    %------------------------------------------------------
                    vbmcBERfullCSIRAvg(iSNR) = ...
                        ((iTrialBlock-1)*vbmcBERfullCSIRAvg(iSNR) + ...
                        mean(vbmcBERfullCSIR))/iTrialBlock;
                    
                    vbmcBERpartCSIRAvg(iSNR) = ...
                        ((iTrialBlock-1)*vbmcBERpartCSIRAvg(iSNR) + ...
                        mean(vbmcBERpartCSIR))/iTrialBlock;
                    
                    vbmcPrSignal_(iTrialBlock) = mean(vbmcPrSignal);
                    vbmcPrNoise_(iTrialBlock) = mean(vbmcPrNoise);
                    vbmcPySignal_(iTrialBlock) = mean(vbmcPySignal);
                    vbmcPyNoise_(iTrialBlock) = mean(vbmcPyNoise);
                    
                    vbmcMaxICIdB_(iTrialBlock) = log10(median(10.^vbmcMaxICIdB));
                    vbmcSIRstatsdB_(iTrialBlock,:) = log10(median(10.^vbmcSIRstatsdB));
                    vbmcSINRstatsdB_(iTrialBlock,:) = log10(median(10.^vbmcSINRstatsdB));
                    vbmcChanJain_(iTrialBlock) = mean(vbmcChanJain);
                    vbmcSirJain_(iTrialBlock) = mean(vbmcSirJain);
                    vbmcSinrJain_(iTrialBlock) = mean(vbmcSinrJain);                    
                    disp(['VBMC: ICImax/SIRmin/avgSIR (dB) = '...
                        num2str(vbmcMaxICIdB_(iTrialBlock)) '/' ...
                        num2str(vbmcSIRstatsdB_(iTrialBlock,1)) '/' ...
                        num2str(vbmcSIRstatsdB_(iTrialBlock,end))]);                                            
                    
                    vbmcPrSignalAvg(iSNR) = mean(vbmcPrSignal_(1:iTrialBlock));
                    vbmcPrNoiseAvg(iSNR) = mean(vbmcPrNoise_(1:iTrialBlock));
                    vbmcPySignalAvg(iSNR) = mean(vbmcPySignal_(1:iTrialBlock));
                    vbmcPyNoiseAvg(iSNR) = mean(vbmcPyNoise_(1:iTrialBlock));
                    vbmcMaxICIdBAvg(iNp,iTau,iAlpha) = ...
                        log10(mean(10.^vbmcMaxICIdB_(1:iTrialBlock)));
                    vbmcSIRstatsdBAvg(iNp,iTau,iAlpha,:) = ...
                        log10(mean(10.^vbmcSIRstatsdB_(1:iTrialBlock,:)));                                        
                    vbmcSINRstatsdBAvg(iNp,iTau,iAlpha,:) = ...
                        log10(mean(10.^vbmcSINRstatsdB_(1:iTrialBlock,:)));                    
                    vbmcSirJainAvg(iNp,iTau,iAlpha) = mean(vbmcSirJain_(1:iTrialBlock));
                    vbmcSinrJainAvg(iNp,iTau,iAlpha) = mean(vbmcSinrJain_(1:iTrialBlock));
                    vbmcChanJainAvg(iNp,iTau,iAlpha) = mean(vbmcChanJain_(1:iTrialBlock));                    
                         
                    %------------------------------------------------------
                    odssBERfullCSIRAvg(iSNR) = ...
                        ((iTrialBlock-1)*odssBERfullCSIRAvg(iSNR) + ...
                        mean(odssBERfullCSIR))/iTrialBlock;
                    
                    odssBERpartCSIRAvg(iSNR) = ...
                        ((iTrialBlock-1)*odssBERpartCSIRAvg(iSNR) + ...
                        mean(odssBERpartCSIR))/iTrialBlock;
                    
                    odssPrSignal_(iTrialBlock) = mean(odssPrSignal);
                    odssPrNoise_(iTrialBlock) = mean(odssPrNoise);
                    odssPySignal_(iTrialBlock) = mean(odssPySignal);
                    odssPyNoise_(iTrialBlock) = mean(odssPyNoise);
                    
                    odssMaxICIdB_(iTrialBlock) = log10(median(10.^odssMaxICIdB));
                    odssSIRstatsdB_(iTrialBlock,:) = log10(median(10.^odssSIRstatsdB));
                    odssSINRstatsdB_(iTrialBlock,:) = log10(median(10.^odssSINRstatsdB));
                    odssChanJain_(iTrialBlock) = mean(odssChanJain);
                    odssSirJain_(iTrialBlock) = mean(odssSirJain);
                    odssSinrJain_(iTrialBlock) = mean(odssSinrJain);                    
                    disp(['ODSS: ICImax/SIRmin/avgSIR (dB) = '...
                        num2str(odssMaxICIdB_(iTrialBlock)) '/' ...
                        num2str(odssSIRstatsdB_(iTrialBlock,1)) '/' ...
                        num2str(odssSIRstatsdB_(iTrialBlock,end))]);                                            
                    
                    odssPrSignalAvg(iSNR) = mean(odssPrSignal_(1:iTrialBlock));
                    odssPrNoiseAvg(iSNR) = mean(odssPrNoise_(1:iTrialBlock));
                    odssPySignalAvg(iSNR) = mean(odssPySignal_(1:iTrialBlock));
                    odssPyNoiseAvg(iSNR) = mean(odssPyNoise_(1:iTrialBlock));
                    odssMaxICIdBAvg(iNp,iTau,iAlpha) = ...
                        log10(mean(10.^odssMaxICIdB_(1:iTrialBlock)));
                    odssSIRstatsdBAvg(iNp,iTau,iAlpha,:) = ...
                        log10(mean(10.^odssSIRstatsdB_(1:iTrialBlock,:)));                                        
                    odssSINRstatsdBAvg(iNp,iTau,iAlpha,:) = ...
                        log10(mean(10.^odssSINRstatsdB_(1:iTrialBlock,:)));                    
                    odssSirJainAvg(iNp,iTau,iAlpha) = mean(odssSirJain_(1:iTrialBlock));
                    odssSinrJainAvg(iNp,iTau,iAlpha) = mean(odssSinrJain_(1:iTrialBlock));
                    odssChanJainAvg(iNp,iTau,iAlpha) = mean(odssChanJain_(1:iTrialBlock));                    
                    
                    %------------------------------------------------------
                    ofdmBERfullCSIRAvg(iSNR) = ...
                        ((iTrialBlock-1)*ofdmBERfullCSIRAvg(iSNR) + ...
                        mean(ofdmBERfullCSIR))/iTrialBlock;
                    
                    ofdmBERpartCSIRAvg(iSNR) = ...
                        ((iTrialBlock-1)*ofdmBERpartCSIRAvg(iSNR) + ...
                        mean(ofdmBERpartCSIR))/iTrialBlock;
                    
                    ofdmPrSignal_(iTrialBlock) = mean(ofdmPrSignal);
                    ofdmPrNoise_(iTrialBlock) = mean(ofdmPrNoise);
                    ofdmPySignal_(iTrialBlock) = mean(ofdmPySignal);
                    ofdmPyNoise_(iTrialBlock) = mean(ofdmPyNoise);        
                    
                    ofdmMaxICIdB_(iTrialBlock) = log10(median(10.^ofdmMaxICIdB));
                    ofdmSIRstatsdB_(iTrialBlock,:) = log10(median(10.^ofdmSIRstatsdB));
                    ofdmSINRstatsdB_(iTrialBlock,:) = log10(median(10.^ofdmSINRstatsdB));
                    ofdmSirJain_(iTrialBlock) = mean(ofdmSirJain);
                    ofdmSinrJain_(iTrialBlock) = mean(ofdmSinrJain);
                    ofdmChanJain_(iTrialBlock) = mean(ofdmChanJain);
                    disp(['OFDM: ICImax/SIRmin/avgSIR (dB) = '...
                        num2str(ofdmMaxICIdB_(iTrialBlock)) '/' ...
                        num2str(ofdmSIRstatsdB_(iTrialBlock,1)) '/' ...
                        num2str(ofdmSIRstatsdB_(iTrialBlock,end))]);                                        
                                        
                    ofdmPrSignalAvg(iSNR) = mean(ofdmPrSignal_(1:iTrialBlock));
                    ofdmPrNoiseAvg(iSNR) = mean(ofdmPrNoise_(1:iTrialBlock));
                    ofdmPySignalAvg(iSNR) = mean(ofdmPySignal_(1:iTrialBlock));
                    ofdmPyNoiseAvg(iSNR) = mean(ofdmPyNoise_(1:iTrialBlock));
                    ofdmMaxICIdBAvg(iNp,iTau,iAlpha) = ...
                        log10(mean(10.^ofdmMaxICIdB_(1:iTrialBlock)));
                    ofdmSIRstatsdBAvg(iNp,iTau,iAlpha,:) = ...
                        log10(mean(10.^ofdmSIRstatsdB_(1:iTrialBlock,:)));                                        
                    ofdmSINRstatsdBAvg(iNp,iTau,iAlpha,:) = ...
                        log10(mean(10.^ofdmSINRstatsdB_(1:iTrialBlock,:)));                    
                    ofdmSirJainAvg(iNp,iTau,iAlpha) = mean(ofdmSirJain_(1:iTrialBlock));
                    ofdmSinrJainAvg(iNp,iTau,iAlpha) = mean(ofdmSinrJain_(1:iTrialBlock));                    
                    ofdmChanJainAvg(iNp,iTau,iAlpha) = mean(ofdmChanJain_(1:iTrialBlock));
                   
                    %------------------------------------------------------
                    otfsBERfullCSIRAvg(iSNR) = ...
                        ((iTrialBlock-1)*otfsBERfullCSIRAvg(iSNR) + ...
                        mean(otfsBERfullCSIR))/iTrialBlock;
                    
                    otfsBERpartCSIRAvg(iSNR) = ...
                        ((iTrialBlock-1)*otfsBERpartCSIRAvg(iSNR) + ...
                        mean(otfsBERpartCSIR))/iTrialBlock;
                    
                    otfsPrSignal_(iTrialBlock) = mean(otfsPrSignal);
                    otfsPrNoise_(iTrialBlock) = mean(otfsPrNoise);
                    otfsPySignal_(iTrialBlock) = mean(otfsPySignal);
                    otfsPyNoise_(iTrialBlock) = mean(otfsPyNoise);                    
                    otfsMaxICIdB_(iTrialBlock) = log10(median(10.^otfsMaxICIdB));
                    otfsSIRstatsdB_(iTrialBlock,:) = log10(median(10.^otfsSIRstatsdB));
                    otfsSINRstatsdB_(iTrialBlock,:) = log10(median(10.^otfsSINRstatsdB));
                    otfsSirJain_(iTrialBlock) = mean(otfsSirJain);
                    otfsSinrJain_(iTrialBlock) = mean(otfsSinrJain);
                    otfsChanJain_(iTrialBlock) = mean(otfsChanJain);     
                    disp(['OTFS: ICImax/SIRmin/avgSIR (dB) = '...
                        num2str(otfsMaxICIdB_(iTrialBlock)) '/' ...
                        num2str(otfsSIRstatsdB_(iTrialBlock,1)) '/' ...
                        num2str(otfsSIRstatsdB_(iTrialBlock,end))]);                                          
                    
                    otfsPrSignalAvg(iSNR) = mean(otfsPrSignal_(1:iTrialBlock));
                    otfsPrNoiseAvg(iSNR) = mean(otfsPrNoise_(1:iTrialBlock));
                    otfsPySignalAvg(iSNR) = mean(otfsPySignal_(1:iTrialBlock));
                    otfsPyNoiseAvg(iSNR) = mean(otfsPyNoise_(1:iTrialBlock));
                    otfsMaxICIdBAvg(iNp,iTau,iAlpha) = ...
                        log10(mean(10.^otfsMaxICIdB_(1:iTrialBlock)));
                    otfsSIRstatsdBAvg(iNp,iTau,iAlpha,:) = ...
                        log10(mean(10.^otfsSIRstatsdB_(1:iTrialBlock,:)));                                        
                    otfsSINRstatsdBAvg(iNp,iTau,iAlpha,:) = ...
                        log10(mean(10.^otfsSINRstatsdB_(1:iTrialBlock,:)));
                    otfsSirJainAvg(iNp,iTau,iAlpha) = mean(otfsSirJain_(1:iTrialBlock));
                    otfsSinrJainAvg(iNp,iTau,iAlpha) = mean(otfsSinrJain_(1:iTrialBlock));                                        
                    otfsChanJainAvg(iNp,iTau,iAlpha) = mean(otfsChanJain_(1:iTrialBlock));
                    %------------------------------------------------------
                    
                    nTrial = iTrialBlock*nTrialsPerBlock;
                    
                    %% Quit Monte-Carlo loop if atleast 100 bit errors are registered for each active decoder
%                     minBitErrs =  ...
%                         min(...
%                         [vbmcBERfullCSIRAvg(iSNR) vbmcBERpartCSIRAvg(iSNR)...
%                         ofdmBERfullCSIRAvg(iSNR) ofdmBERpartCSIRAvg(iSNR)...
%                         otfsBERfullCSIRAvg(iSNR) otfsBERpartCSIRAvg(iSNR)]*...
%                         nTrial*N ...
%                         );

                    minBitErrs =  ...
                        min(...
                        [vbmcBERfullCSIRAvg(iSNR) vbmcBERpartCSIRAvg(iSNR)...
                        odssBERfullCSIRAvg(iSNR) odssBERpartCSIRAvg(iSNR)...
                        ofdmBERfullCSIRAvg(iSNR) ofdmBERpartCSIRAvg(iSNR)...
                        otfsBERpartCSIRAvg(iSNR) otfsBERfullCSIRAvg(iSNR)]*... %otfsBERfullCSIRAvg(iSNR): avoided to speed up trials!!!
                        nTrial*N ...
                        );

                    if ((minBitErrs > NBITERRMIN) && (nTrial > NtrialMin)) || (toc(tDispl) > displayRate*60)
                        disp(' ');
                        disp(['VBMC Trial # ' num2str(nTrial) ': ' num2str([vbmcBERfullCSIRAvg(iSNR) vbmcBERpartCSIRAvg(iSNR)])...
                            ' (SNR = ' num2str(SNRlist(iSNR)) ' dB)' ' [bitErrs = ' num2str(minBitErrs) ']']);
                        disp(['ODSS Trial # ' num2str(nTrial) ': ' num2str([odssBERfullCSIRAvg(iSNR) odssBERpartCSIRAvg(iSNR)])...
                            ' (SNR = ' num2str(SNRlist(iSNR)) ' dB)' ' [bitErrs = ' num2str(minBitErrs) ']']);                        
                        disp(['OFDM Trial # ' num2str(nTrial) ': ' num2str([ofdmBERfullCSIRAvg(iSNR) ofdmBERpartCSIRAvg(iSNR)])...
                            ' (SNR = ' num2str(SNRlist(iSNR)) ' dB)' ' [bitErrs = ' num2str(minBitErrs) ']']);
                        disp(['OTFS Trial # ' num2str(nTrial) ': ' num2str([otfsBERfullCSIRAvg(iSNR) otfsBERpartCSIRAvg(iSNR)])...
                            ' (SNR = ' num2str(SNRlist(iSNR)) ' dB)' ' [bitErrs = ' num2str(minBitErrs) ']']);
                        
                        if ((minBitErrs > NBITERRMIN) && (nTrial > NtrialMin))
                            disp(['Monte-Carlo loop terminated since total bit error count > '...
                                num2str(NBITERRMIN) ' for all active decoders and # of trials='...
                                num2str(nTrial) ' > ' num2str(NtrialMin)]);
                        else
                            tDispl  = tic; %reset the timer for displaying results
                        end
                        
                        % Display Jain indices for each scheme...
                        disp([' alphaMax = ' num2str(alphaMax) ','...
                            ' tauMax = ' num2str(tauMax*1e3) ' ms,'...
                            ' Np = ' num2str(Np) ' paths']);                        
                        disp(['VBMC: Mean Jain Index = ' num2str(vbmcChanJainAvg(iNp,iTau,iAlpha))]);
                        disp(['ODSS: Mean Jain Index = ' num2str(odssChanJainAvg(iNp,iTau,iAlpha))]);
                        disp(['OFDM: Mean Jain Index = ' num2str(ofdmChanJainAvg(iNp,iTau,iAlpha))]);
                        disp(['OTFS: Mean Jain Index = ' num2str(otfsChanJainAvg(iNp,iTau,iAlpha))]);
                        disp(['(iNp,iTau,iAlpha) : (' num2str(iNp) ', ' ...
                            num2str(iTau) ', ' num2str(iAlpha) ')']);
                        %---------------------------------------------------
                        %Update plots before proceeding to next SNR / every displayRate minutes
                        % F.1) Plot bit error rates
                        % berFig = 18
                        set(0, 'CurrentFigure', berFig) , cla,
                        semilogy(SNRlist(1:iSNR), vbmcBERfullCSIRAvg(1:iSNR), 'b-*'), hold on, grid on
                        semilogy(SNRlist(1:iSNR), vbmcBERpartCSIRAvg(1:iSNR), 'r-o'), hold on, grid on
                        semilogy(SNRlist(1:iSNR), odssBERfullCSIRAvg(1:iSNR), 'g-<'), hold on, grid on
                        semilogy(SNRlist(1:iSNR), odssBERpartCSIRAvg(1:iSNR), 'm-.d'), hold on, grid on                        
                        semilogy(SNRlist(1:iSNR), ofdmBERfullCSIRAvg(1:iSNR), 'b-.x'), hold on, grid on
                        semilogy(SNRlist(1:iSNR), ofdmBERpartCSIRAvg(1:iSNR), 'r-.s'), hold on, grid on
                        semilogy(SNRlist(1:iSNR), otfsBERfullCSIRAvg(1:iSNR), 'b-.>'), hold on, grid on
                        semilogy(SNRlist(1:iSNR), otfsBERpartCSIRAvg(1:iSNR), 'r-.p'), hold on, grid on                        
                        xlabel('SNR (dB)'), ylabel('BER'),
                        legend('VBMC (Full CSIR)', 'VBMC (Partial CSIR)',...
                            'ODSS (Full CSIR)', 'ODSS (Partial CSIR)',...
                            'OFDM (Full CSIR)', 'OFDM (Partial CSIR)',...
                            'OTFS (Full CSIR)', 'OTFS (Partial CSIR)')
                        title([num2str(nTrial) '/' num2str(Ntrials) ' trials, bit errors (best scheme) = ' num2str(minBitErrs)])
                        
                        % Plot measured VBMC signal, symbol, noise powers and SNR
                        % vbmcMeasFig = 99
                        set(0, 'CurrentFigure', vbmcMeasFig) , cla, %vbmcMeasFig = 99;
                        subplot(221),
                        plot(1:iTrialBlock, 10*log10(vbmcPrSignal_(1:iTrialBlock)), '-*'), hold on, grid on
                        plot(1:iTrialBlock, ones(iTrialBlock,1)*10*log10(Ps0V), 'k'), hold off %overlay expected signal power
                        xlabel('trial block #'), ylabel('10 log_{10} P_{signal}'), legend('Estimated', 'Expected')
                        title(['VBMC: Rx Wav Pow (dB) = '...
                            num2str(round(10*log10(vbmcPrSignalAvg(iSNR))*10)/10) ' ('...
                            num2str(round(10*log10(Ps0V)*10)/10) ')'])
                        subplot(222),
                        plot(1:iTrialBlock, 10*log10(vbmcPrSignal_(1:iTrialBlock)./vbmcPrNoise_(1:iTrialBlock))+10*log10(FsV/B), '-s'),
                        hold on, grid on, plot(1:iTrialBlock, ones(iTrialBlock,1)*(SNRv+10*log10(FsV/B)), 'k'), hold off %overlay expected SNR
                        xlabel('trial block #'), ylabel('SNR (dB)'), legend('Estimated', 'Expected')
                        title(['VBMC: Rx Wav SNR (dB) = '...
                            num2str(round((10*log10(vbmcPrSignalAvg(iSNR)./vbmcPrNoiseAvg(iSNR))+10*log10(FsV/B))*10)/10)...
                            ' (' num2str(round((SNRv+10*log10(FsV/B))*10)/10) ')'])
                        subplot(223),
                        plot(1:iTrialBlock, 10*log10(vbmcPyNoise_(1:iTrialBlock)/N), '-*'), hold on, grid on
                        plot(1:iTrialBlock, ones(iTrialBlock,1)*10*log10(mean(PnoisetV)), 'k'), hold off %overlay expected signal power
                        xlabel('trial block #'), ylabel('10 log_{10} P_{noise,symb}'), legend('Estimated', 'Expected')
                        title(['VBMC: Rx Symb Noise Pow (dB) = '...
                            num2str(round(10*log10(vbmcPyNoiseAvg(iSNR)/N)*10)/10) ' ('...
                            num2str(round(10*log10(mean(PnoisetV))*10)/10) ')'])
                        subplot(224),
                        plot(1:iTrialBlock, 10*log10(vbmcPySignal_(1:iTrialBlock)./vbmcPyNoise_(1:iTrialBlock)), '-s'), grid on
                        %             hold on, plot(1:trial, ones(trial,1)*SNR, 'k'), hold off %overlay expected SNR
                        xlabel('trial block #'), ylabel('SNR (dB)'), legend('Estimated')%, 'Expected')
                        title(['VBMC: Rx Symb SNR (dB) = '...
                            num2str(round(10*log10(vbmcPySignalAvg(iSNR)./vbmcPyNoiseAvg(iSNR))*10)/10)]);%...
                        %                 ' (' num2str(round(SNR*10)/10) ')'])
                        
                        % Plot measured ODSS signal, symbol, noise powers and SNR
                        % odssMeasFig = 99999
                        set(0, 'CurrentFigure', odssMeasFig) , cla, %odssMeasFig = 99999;
                        subplot(221),
                        plot(1:iTrialBlock, 10*log10(odssPrSignal_(1:iTrialBlock)), '-*'), hold on, grid on
                        plot(1:iTrialBlock, ones(iTrialBlock,1)*10*log10(Ps0V), 'k'), hold off %overlay expected signal power
                        xlabel('trial block #'), ylabel('10 log_{10} P_{signal}'), legend('Estimated', 'Expected')
                        title(['ODSS: Rx Wav Pow (dB) = '...
                            num2str(round(10*log10(odssPrSignalAvg(iSNR))*10)/10) ' ('...
                            num2str(round(10*log10(Ps0V)*10)/10) ')'])
                        subplot(222),
                        plot(1:iTrialBlock, 10*log10(odssPrSignal_(1:iTrialBlock)./odssPrNoise_(1:iTrialBlock))+10*log10(FsV/B), '-s'),
                        hold on, grid on, plot(1:iTrialBlock, ones(iTrialBlock,1)*(SNRv+10*log10(FsV/B)), 'k'), hold off %overlay expected SNR
                        xlabel('trial block #'), ylabel('SNR (dB)'), legend('Estimated', 'Expected')
                        title(['ODSS: Rx Wav SNR (dB) = '...
                            num2str(round((10*log10(odssPrSignalAvg(iSNR)./odssPrNoiseAvg(iSNR))+10*log10(FsV/B))*10)/10)...
                            ' (' num2str(round((SNRv+10*log10(FsV/B))*10)/10) ')'])
                        subplot(223),
                        plot(1:iTrialBlock, 10*log10(odssPyNoise_(1:iTrialBlock)/N), '-*'), hold on, grid on
                        plot(1:iTrialBlock, ones(iTrialBlock,1)*10*log10(mean(PnoisetD)), 'k'), hold off %overlay expected signal power
                        xlabel('trial block #'), ylabel('10 log_{10} P_{noise,symb}'), legend('Estimated', 'Expected')
                        title(['ODSS: Rx Symb Noise Pow (dB) = '...
                            num2str(round(10*log10(odssPyNoiseAvg(iSNR)/N)*10)/10) ' ('...
                            num2str(round(10*log10(mean(PnoisetD))*10)/10) ')'])
                        subplot(224),
                        plot(1:iTrialBlock, 10*log10(odssPySignal_(1:iTrialBlock)./odssPyNoise_(1:iTrialBlock)), '-s'), grid on
                        %             hold on, plot(1:trial, ones(trial,1)*SNR, 'k'), hold off %overlay expected SNR
                        xlabel('trial block #'), ylabel('SNR (dB)'), legend('Estimated')%, 'Expected')
                        title(['ODSS: Rx Symb SNR (dB) = '...
                            num2str(round(10*log10(odssPySignalAvg(iSNR)./odssPyNoiseAvg(iSNR))*10)/10)]);%...
                        %                 ' (' num2str(round(SNR*10)/10) ')'])                        
                        
                        % Plot measured OFDM signal, symbol, noise powers and SNR
                        % ofdmMeasFig = 999
                        set(0, 'CurrentFigure', ofdmMeasFig) , cla, % ofdmMeasFig = 999;
                        subplot(221),
                        plot(1:iTrialBlock, 10*log10(ofdmPrSignal_(1:iTrialBlock)), '-*'), hold on, grid on
                        plot(1:iTrialBlock, ones(iTrialBlock,1)*10*log10(Ps0O), 'k'), hold off %overlay expected signal power
                        xlabel('trial block #'), ylabel('10 log_{10} P_{signal}'), legend('Estimated', 'Expected')
                        title(['OFDM: Rx Wav Pow (dB) = '...
                            num2str(round(10*log10(mean(ofdmPrSignalAvg(iSNR)))*10)/10) ' ('...
                            num2str(round(10*log10(Ps0O)*10)/10) ')'])
                        subplot(222),
                        plot(1:iTrialBlock, 10*log10(ofdmPrSignal_(1:iTrialBlock)./ofdmPrNoise_(1:iTrialBlock))+10*log10(FsO/B), '-s'),
                        hold on, grid on, plot(1:iTrialBlock, ones(iTrialBlock,1)*(SNRo+10*log10(FsO/B)), 'k'), hold off %overlay expected SNR
                        xlabel('trial block #'), ylabel('SNR (dB)'), legend('Estimated', 'Expected')
                        title(['OFDM: Rx Wav SNR (dB) = '...
                            num2str(round((10*log10(ofdmPrSignalAvg(iSNR)./ofdmPrNoiseAvg(iSNR))+10*log10(FsO/B))*10)/10)...
                            ' (' num2str(round((SNRo+10*log10(FsO/B))*10)/10) ')'])
                        subplot(223),
                        plot(1:iTrialBlock, 10*log10(ofdmPyNoise_(1:iTrialBlock)/N), '-*'), hold on, grid on
                        plot(1:iTrialBlock, ones(iTrialBlock,1)*10*log10(mean(PnoisetO)), 'k'), hold off %overlay expected signal power
                        xlabel('trial block #'), ylabel('10 log_{10} P_{noise,symb}'), legend('Estimated', 'Expected')
                        title(['OFDM: Rx Symb Noise Pow (dB) = '...
                            num2str(round(10*log10(ofdmPyNoiseAvg(iSNR)/N)*10)/10) ' ('...
                            num2str(round(10*log10(mean(PnoisetO))*10)/10) ')'])
                        subplot(224),
                        plot(1:iTrialBlock, 10*log10(ofdmPySignal_(1:iTrialBlock)./ofdmPyNoise_(1:iTrialBlock)), '-s'), grid on
                        %             hold on, plot(1:trial, ones(trial,1)*SNR, 'k'), hold off %overlay expected SNR
                        xlabel('trial block #'), ylabel('SNR (dB)'), legend('Estimated')%, 'Expected')
                        title(['OFDM: Rx Symb SNR (dB) = '...
                            num2str(round(10*log10(ofdmPySignalAvg(iSNR)./ofdmPyNoiseAvg(iSNR))*10)/10)]);%...
                        %                 ' (' num2str(round(SNR*10)/10) ')'])
                        
                        % Plot measured OTFS signal, symbol, noise powers and SNR
                        % otfsMeasFig = 9999
                        set(0, 'CurrentFigure', otfsMeasFig) , cla, % otfsMeasFig = 999;
                        subplot(221),
                        plot(1:iTrialBlock, 10*log10(otfsPrSignal_(1:iTrialBlock)), '-*'), hold on, grid on
                        plot(1:iTrialBlock, ones(iTrialBlock,1)*10*log10(Ps0O), 'k'), hold off %overlay expected signal power
                        xlabel('trial block #'), ylabel('10 log_{10} P_{signal}'), legend('Estimated', 'Expected')
                        title(['OTFS: Rx Wav Pow (dB) = '...
                            num2str(round(10*log10(mean(otfsPrSignal_(1:iTrialBlock)))*10)/10) ' ('...
                            num2str(round(10*log10(Ps0O)*10)/10) ')'])
                        subplot(222),
                        plot(1:iTrialBlock, 10*log10(otfsPrSignal_(1:iTrialBlock)./otfsPrNoise_(1:iTrialBlock))+10*log10(FsO/B), '-s'),
                        hold on, grid on, plot(1:iTrialBlock, ones(iTrialBlock,1)*(SNRo+10*log10(FsO/B)), 'k'), hold off %overlay expected SNR
                        xlabel('trial block #'), ylabel('SNR (dB)'), legend('Estimated', 'Expected')
                        title(['OTFS: Rx Wav SNR (dB) = '...
                            num2str(round((10*log10(otfsPrSignalAvg(iSNR)./otfsPrNoiseAvg(iSNR))+10*log10(FsO/B))*10)/10)...
                            ' (' num2str(round((SNRo+10*log10(FsO/B))*10)/10) ')'])
                        subplot(223),
                        plot(1:iTrialBlock, 10*log10(otfsPyNoise_(1:iTrialBlock)/N), '-*'), hold on, grid on
                        plot(1:iTrialBlock, ones(iTrialBlock,1)*10*log10(mean(PnoisetO)), 'k'), hold off %overlay expected signal power
                        xlabel('trial block #'), ylabel('10 log_{10} P_{noise,symb}'), legend('Estimated', 'Expected')
                        title(['OTFS: Rx Symb Noise Pow (dB) = '...
                            num2str(round(10*log10(otfsPyNoiseAvg(iSNR)/N)*10)/10) ' ('...
                            num2str(round(10*log10(mean(PnoisetO))*10)/10) ')'])
                        subplot(224),
                        plot(1:iTrialBlock, 10*log10(otfsPySignal_(1:iTrialBlock)./otfsPyNoise_(1:iTrialBlock)), '-s'), grid on
                        %             hold on, plot(1:trial, ones(trial,1)*SNR, 'k'), hold off %overlay expected SNR
                        xlabel('trial block #'), ylabel('SNR (dB)'), legend('Estimated')%, 'Expected')
                        title(['OTFS: Rx Symb SNR (dB) = '...
                            num2str(round(10*log10(otfsPySignalAvg(iSNR)./otfsPyNoiseAvg(iSNR))*10)/10)]);%...
                        %                 ' (' num2str(round(SNR*10)/10) ')'])
                        
                        % F.2a) Plot measured VBMC Rx Channel Matrix
                        % vbmcChanMatFig = 100
                        set(0, 'CurrentFigure', vbmcChanMatFig) , cla, % vbmcChanMatFig = 100;
                        imagesc(10*log10(abs(HtV(:,:,1)))), colorbar
                        xlabel('subcarrier index (n)'), ylabel('subcarrier index (m)'),
                        title(['VBMC Rx: ICI_{max} = '...
                            num2str(round(log10(10.^vbmcMaxICIdBAvg(iNp,iTau,iAlpha))*10)/10)...
                            ' dB, SIR_{min} = '...
                            num2str(round(log10(10.^vbmcSIRstatsdBAvg(iNp,iTau,iAlpha,1))*10)/10) ' dB']);                        
                        
                        % F.2b) Plot measured ODSS Rx Channel Matrix
                        % odssChanMatFig = 100000
                        set(0, 'CurrentFigure', odssChanMatFig) , cla, % odssChanMatFig = 10000;
                        imagesc(10*log10(abs(HtV(:,:,1)))), colorbar
                        xlabel('subcarrier index (n)'), ylabel('subcarrier index (m)'),
                        title(['ODSS Rx: ICI_{max} = '...
                            num2str(round(log10(10.^odssMaxICIdBAvg(iNp,iTau,iAlpha))*10)/10)...
                            ' dB, SIR_{min} = '...
                            num2str(round(log10(10.^odssSIRstatsdBAvg(iNp,iTau,iAlpha,1))*10)/10) ' dB']);                        
                        
                        % F.3a) Plot measured OFDM Rx Channel Matrix
                        % ofdmChanMatFig = 1000
                        set(0, 'CurrentFigure', ofdmChanMatFig) , cla, % ofdmChanMatFig = 1000;
                        imagesc(10*log10(abs(HtO(:,:,1)))), colorbar
                        xlabel('subcarrier index (n)'), ylabel('subcarrier index (m)'),
                        title(['OFDM Rx: ICI_{max} = '...
                            num2str(round(log10(10.^ofdmMaxICIdBAvg(iNp,iTau,iAlpha))*10)/10)...
                            ' dB, SIR_{min} = '...
                            num2str(round(log10(10.^ofdmSIRstatsdBAvg(iNp,iTau,iAlpha,1))*10)/10) ' dB']);
                        
                        % F.3b) Plot measured OTFS Rx Channel Matrix
                        % otfsChanMatFig = 10000
                        set(0, 'CurrentFigure', otfsChanMatFig) , cla, % otfsChanMatFig = 10000;
                        imagesc(10*log10(abs(HtT(:,:,1)))), colorbar
                        xlabel('subcarrier index (n)'), ylabel('subcarrier index (m)'),
                        title(['OTFS Rx: ICI_{max} = '...
                            num2str(round(log10(10.^otfsMaxICIdBAvg(iNp,iTau,iAlpha))*10)/10)...
                            ' dB, SIR_{min} = '...
                            num2str(round(log10(10.^otfsSIRstatsdBAvg(iNp,iTau,iAlpha,1))*10)/10) ' dB']);                        
                        
                        pause(0.01);
                        
                        if toc(tSave) > saveRate*60 %save results on elapse of saveRate minutes
                            %% Save results & plots before the run for next alphaMax
                            simParamStr = ['_tauMax' num2str(tauMax*1e3) 'ms'...
                                '_alphaMax' num2str(alphaMax)...
                                '_alphaMaxDes' num2str(alphaMaxDes) ...
                                '_Np' num2str(Np)...
                                '_T' num2str(T*1e3) 'ms'...
                                '_N' num2str(N)...
                                '_SNR_' num2str(SNRlist(1)) 'To' num2str(SNRlist(end)) 'dB' ...
                                '_NTRIALS' num2str(Ntrials)...
                                '_NBITERRMIN' num2str(NBITERRMIN)];
                            
                            save(['figures/vbmcVofdmVotfsVodss_simResults' simParamStr '.mat'], ...
                                'NpList', 'tauMaxList', 'alphaMaxList', 'SNRlist',...
                                'tauMax', 'alphaMax', 'Np',...
                                'fL', 'fH', 'dF', 'betaB', 'T', 'B', 'fc', ...
                                'alphaMaxDes', 'N', 'FsV',...
                                'Nfft', 'FsO', 'Lcp',...
                                'betaModn',...
                                'Ntrials', 'NBITERRMIN', 'NtrialMin',...
                                'vbmcBERfullCSIRAvg', 'vbmcBERpartCSIRAvg',...
                                'vbmcPrSignalAvg', 'vbmcPrSigCalAvg', 'vbmcPySignalAvg',...
                                'vbmcPrNoiseAvg', 'vbmcPyNoiseAvg',...
                                'vbmcMaxICIdBAvg', 'vbmcSIRstatsdBAvg',...
                                'vbmcSirJainAvg', 'vbmcSinrJainAvg', 'vbmcChanJainAvg',...
                                'odssBERfullCSIRAvg', 'odssBERpartCSIRAvg',...
                                'odssPrSignalAvg', 'odssPrSigCalAvg', 'odssPySignalAvg',...
                                'odssPrNoiseAvg', 'odssPyNoiseAvg',...
                                'odssMaxICIdBAvg', 'odssSIRstatsdBAvg',...
                                'odssSirJainAvg', 'odssSinrJainAvg', 'odssChanJainAvg',...                                
                                'ofdmBERfullCSIRAvg', 'ofdmBERpartCSIRAvg',...
                                'ofdmPrSignalAvg', 'ofdmPrSigCalAvg', 'ofdmPySignalAvg',...
                                'ofdmPrNoiseAvg', 'ofdmPyNoiseAvg',...
                                'ofdmMaxICIdBAvg', 'ofdmSIRstatsdBAvg',...
                                'ofdmSirJainAvg', 'ofdmSinrJainAvg', 'ofdmChanJainAvg',...
                                'otfsBERfullCSIRAvg', 'otfsBERpartCSIRAvg',...
                                'otfsPrSignalAvg', 'otfsPrSigCalAvg', 'otfsPySignalAvg',...
                                'otfsPrNoiseAvg', 'otfsPyNoiseAvg',...
                                'otfsMaxICIdBAvg', 'otfsSIRstatsdBAvg',...
                                'otfsSirJainAvg', 'otfsSinrJainAvg', 'otfsChanJainAvg');                                
                            
                            set(0, 'CurrentFigure', berFig), %figure(18),
                            savefig(['figures/vbmcVofdmVotfsVodss_berVsnr' simParamStr '.fig']);
                            pause(0.1)
                            
                            set(0, 'CurrentFigure', vbmcChanMatFig), %figure(100),
                            savefig(['figures/vbmc_RxChanMat' simParamStr '.fig']);
                            pause(0.1)
                            
                            set(0, 'CurrentFigure', ofdmChanMatFig), %figure(1000),
                            savefig(['figures/ofdm_RxChanMat' simParamStr '.fig']);
                            pause(0.1)
                            
                            set(0, 'CurrentFigure', otfsChanMatFig), %figure(10000),
                            savefig(['figures/otfs_RxChanMat' simParamStr '.fig']);
                            pause(0.1)                            
                            
                            set(0, 'CurrentFigure', symbSnrFig), %figure(27),
                            savefig(['figures/vbmcVofdmVotfsVodss_SymbSNR' simParamStr '.fig']);
                            pause(0.1)
                            
                            tSave  = tic; %reset the timer for saving results
                        end
                        
                        if ((minBitErrs > NBITERRMIN) && (nTrial > NtrialMin))
                            break;
                        end
                    end
                    
                end
                
                %----------------------------------------------------------------------
                vbmcPrSignalAvg(iSNR) = mean(vbmcPrSignal_);
                vbmcPrNoiseAvg(iSNR) = mean(vbmcPrNoise_);
                vbmcPySignalAvg(iSNR) = mean(vbmcPySignal_);
                vbmcPyNoiseAvg(iSNR) = mean(vbmcPyNoise_);
%                 vbmcMaxICIdBAvg(iNp,iTau,iAlpha) = ...
%                     log10(mean(10.^vbmcMaxICIdB_(1:iTrialBlock)));
                %----------------------------------------------------------------------
                odssPrSignalAvg(iSNR) = mean(odssPrSignal_);
                odssPrNoiseAvg(iSNR) = mean(odssPrNoise_);
                odssPySignalAvg(iSNR) = mean(odssPySignal_);
                odssPyNoiseAvg(iSNR) = mean(odssPyNoise_);
%                 odssMaxICIdBAvg(iNp,iTau,iAlpha) = ...
%                     log10(mean(10.^odssMaxICIdB_(1:iTrialBlock)));
                %----------------------------------------------------------
                ofdmPrSignalAvg(iSNR) = mean(ofdmPrSignal_);
                ofdmPrNoiseAvg(iSNR) = mean(ofdmPrNoise_);
                ofdmPySignalAvg(iSNR) = mean(ofdmPySignal_);
                ofdmPyNoiseAvg(iSNR) = mean(ofdmPyNoise_);
%                 ofdmMaxICIdBAvg(iNp,iTau,iAlpha) = ...
%                     log10(mean(10.^ofdmMaxICIdB_(1:iTrialBlock)));
                %----------------------------------------------------------
                otfsPrSignalAvg(iSNR) = mean(otfsPrSignal_);
                otfsPrNoiseAvg(iSNR) = mean(otfsPrNoise_);
                otfsPySignalAvg(iSNR) = mean(otfsPySignal_);
                otfsPyNoiseAvg(iSNR) = mean(otfsPyNoise_);
%                 otfsMaxICIdBAvg(iNp,iTau,iAlpha) = ...
%                     log10(mean(10.^otfsMaxICIdB_(1:iTrialBlock)));
                %----------------------------------------------------------
                
                % Update & Display Jain indices for each scheme...                
                disp([' alphaMax = ' num2str(alphaMax) ','...
                    ' tauMax = ' num2str(tauMax*1e3) ' ms,'...
                    ' Np = ' num2str(Np) ' paths']);  
                vbmcChanJainAvg(iNp,iTau,iAlpha) = mean(vbmcChanJain_(1:iTrialBlock));
                disp(['VBMC: Mean Jain Index = ' num2str(vbmcChanJainAvg(iNp,iTau,iAlpha))]);
                odssChanJainAvg(iNp,iTau,iAlpha) = mean(odssChanJain_(1:iTrialBlock));
                disp(['ODSS: Mean Jain Index = ' num2str(odssChanJainAvg(iNp,iTau,iAlpha))]);                
                ofdmChanJainAvg(iNp,iTau,iAlpha) = mean(ofdmChanJain_(1:iTrialBlock));
                disp(['OFDM: Mean Jain Index = ' num2str(ofdmChanJainAvg(iNp,iTau,iAlpha))]);
                otfsChanJainAvg(iNp,iTau,iAlpha) = mean(otfsChanJain_(1:iTrialBlock));
                disp(['OTFS: Mean Jain Index = ' num2str(otfsChanJainAvg(iNp,iTau,iAlpha))]);                           
                disp(['(iNp,iTau,iAlpha) : (' num2str(iNp) ', ' ...
                    num2str(iTau) ', ' num2str(iAlpha) ')']);                
                %----------------------------------------------------------
                % F.1) Plot bit error rates
                % berFig = 18
                set(0, 'CurrentFigure', berFig) , cla, % berFig = 18
                semilogy(SNRlist(1:iSNR), vbmcBERfullCSIRAvg(1:iSNR), 'b-*'), hold on, grid on
                semilogy(SNRlist(1:iSNR), vbmcBERpartCSIRAvg(1:iSNR), 'r-o'), hold on, grid on
                semilogy(SNRlist(1:iSNR), odssBERfullCSIRAvg(1:iSNR), 'g-<'), hold on, grid on
                semilogy(SNRlist(1:iSNR), odssBERpartCSIRAvg(1:iSNR), 'm-.'), hold on, grid on                
                semilogy(SNRlist(1:iSNR), ofdmBERfullCSIRAvg(1:iSNR), 'b-.x'), hold on, grid on
                semilogy(SNRlist(1:iSNR), ofdmBERpartCSIRAvg(1:iSNR), 'r-.s'), hold on, grid on
                semilogy(SNRlist(1:iSNR), otfsBERfullCSIRAvg(1:iSNR), 'b-.>'), hold on, grid on
                semilogy(SNRlist(1:iSNR), otfsBERpartCSIRAvg(1:iSNR), 'r-.p'), hold on, grid on
                xlabel('SNR (dB)'), ylabel('BER'),
                legend('VBMC (Full CSIR)', 'VBMC (Partial CSIR)',...
                    'ODSS (Full CSIR)', 'ODSS (Partial CSIR)',...
                    'OFDM (Full CSIR)', 'OFDM (Partial CSIR)',...
                    'OTFS (Full CSIR)', 'OTFS (Partial CSIR)')
                title(['\tau_{max} = ' num2str(tauMax*1e3)...
                    ' ms, \alpha_{max} = '  num2str(alphaMax)...
                    ', Np = ' num2str(Np)]);
                
                % F.4) Plot symbol SNR vs wav. SNR at the Rx for VBMC/OFDM
                % symbSnrFig = 27
                set(0, 'CurrentFigure', symbSnrFig), cla, % symbSnrFig = 27
                plot(SNRlist(1:iSNR),...
                    10*log10(vbmcPySignalAvg(1:iSNR)./vbmcPyNoiseAvg(1:iSNR)), '-o'),
                hold on, grid on,
                plot(SNRlist(1:iSNR),...
                    10*log10(odssPySignalAvg(1:iSNR)./odssPyNoiseAvg(1:iSNR)), '-d'),                
                plot(SNRlist(1:iSNR),...
                    10*log10(ofdmPySignalAvg(1:iSNR)./ofdmPyNoiseAvg(1:iSNR)), '-s'),
                plot(SNRlist(1:iSNR),...
                    10*log10(otfsPySignalAvg(1:iSNR)./otfsPyNoiseAvg(1:iSNR)), '->'),                
                xlabel('Rx Waveform SNR (dB)'), ylabel('Rx Symbol SNR (dB)'),
                legend('VBMC','ODSS','OFDM','OTFS')
                title( ['\tau_{max} = ' num2str(tauMax*1e3)...
                    ' ms, \alpha_{max} = '  num2str(alphaMax)...
                    ', Np = ' num2str(Np) ...
                    ': ICI_{max} VBMC/ODSS/OFDM/OTFS = '...
                    num2str(round(log10(mean(10.^vbmcMaxICIdBAvg(iNp,iTau,iAlpha)))*10)/10)...
                    ' / '...
                    num2str(round(log10(mean(10.^odssMaxICIdBAvg(iNp,iTau,iAlpha)))*10)/10)...
                    ' / '...                    
                    num2str(round(log10(mean(10.^ofdmMaxICIdBAvg(iNp,iTau,iAlpha)))*10)/10)...
                    ' / '...
                    num2str(round(log10(mean(10.^otfsMaxICIdBAvg(iNp,iTau,iAlpha)))*10)/10)] );
                
                pause(0.01); %wait a bit for all updates to finish...
                
                disp(['%%%%%%%%% Run Time = ' num2str(toc(tExecSnr)) 's:'...
                    ' SNR = ' num2str(SNRlist(iSNR)) ' dB,'...
                    ' alphaMax = ' num2str(alphaMax) ','...
                    ' tauMax = ' num2str(tauMax*1e3) ' ms,'...
                    ' Np = ' num2str(Np) ' paths %%%%%%%%%']);
            end %end of SNR loop
            %----------------------------------------------------------------------
            
            %% Save results & plots before the run for next alphaMax
            simParamStr = ['_tauMax' num2str(tauMax*1e3) 'ms'...
                '_alphaMax' num2str(alphaMax)...
                '_alphaMaxDes' num2str(alphaMaxDes) ...
                '_Np' num2str(Np)...
                '_T' num2str(T*1e3) 'ms'...
                '_N' num2str(N)...
                '_SNR_' num2str(SNRlist(1)) 'To' num2str(SNRlist(end)) 'dB' ...
                '_NTRIALS' num2str(Ntrials)...
                '_NBITERRMIN' num2str(NBITERRMIN)];
            
            save(['figures/vbmcVofdmVotfsVodss_simResults' simParamStr '.mat'], ...
                'NpList', 'tauMaxList', 'alphaMaxList', 'SNRlist',...
                'tauMax', 'alphaMax', 'Np',...
                'fL', 'fH', 'dF', 'betaB', 'T', 'B', 'fc', ...
                'alphaMaxDes', 'N', 'FsV',...
                'Nfft', 'FsO', 'Lcp',...
                'betaModn',...
                'Ntrials', 'NBITERRMIN', 'NtrialMin',...
                'vbmcBERfullCSIRAvg', 'vbmcBERpartCSIRAvg',...
                'vbmcPrSignalAvg', 'vbmcPrSigCalAvg', 'vbmcPySignalAvg',...
                'vbmcPrNoiseAvg', 'vbmcPyNoiseAvg',...
                'vbmcMaxICIdBAvg', 'vbmcSIRstatsdBAvg', 'vbmcSINRstatsdBAvg',...
                'vbmcSirJainAvg', 'vbmcSinrJainAvg', 'vbmcChanJainAvg',...
                'odssBERfullCSIRAvg', 'odssBERpartCSIRAvg',...
                'odssPrSignalAvg', 'odssPrSigCalAvg', 'odssPySignalAvg',...
                'odssPrNoiseAvg', 'odssPyNoiseAvg',...
                'odssMaxICIdBAvg', 'odssSIRstatsdBAvg', 'odssSINRstatsdBAvg',...
                'odssSirJainAvg', 'odssSinrJainAvg', 'odssChanJainAvg',...                
                'ofdmBERfullCSIRAvg', 'ofdmBERpartCSIRAvg',...
                'ofdmPrSignalAvg', 'ofdmPrSigCalAvg', 'ofdmPySignalAvg',...
                'ofdmPrNoiseAvg', 'ofdmPyNoiseAvg',...
                'ofdmMaxICIdBAvg', 'ofdmSIRstatsdBAvg', 'ofdmSINRstatsdBAvg',...
                'ofdmSirJainAvg', 'ofdmSinrJainAvg', 'ofdmChanJainAvg',...
                'otfsBERfullCSIRAvg', 'otfsBERpartCSIRAvg',...
                'otfsPrSignalAvg', 'otfsPrSigCalAvg', 'otfsPySignalAvg',...
                'otfsPrNoiseAvg', 'otfsPyNoiseAvg',...
                'otfsMaxICIdBAvg', 'otfsSIRstatsdBAvg', 'otfsSINRstatsdBAvg',...
                'otfsSirJainAvg', 'otfsSinrJainAvg', 'otfsChanJainAvg');
            
            set(0, 'CurrentFigure', berFig), %figure(18),
            savefig(['figures/vbmcVofdmVotfsVodss_berVsnr' simParamStr '.fig']);
            pause(0.1)
            
            set(0, 'CurrentFigure', vbmcChanMatFig), %figure(100),
            savefig(['figures/vbmc_RxChanMat' simParamStr '.fig']);
            pause(0.1)
            
            set(0, 'CurrentFigure', odssChanMatFig), %figure(100),
            savefig(['figures/vbmc_RxChanMat' simParamStr '.fig']);
            pause(0.1)            
            
            set(0, 'CurrentFigure', ofdmChanMatFig), %figure(1000),
            savefig(['figures/ofdm_RxChanMat' simParamStr '.fig']);
            pause(0.1)
            
            set(0, 'CurrentFigure', otfsChanMatFig), %figure(10000),
            savefig(['figures/otfs_RxChanMat' simParamStr '.fig']);
            pause(0.1)
            
            set(0, 'CurrentFigure', symbSnrFig), %figure(27),
            savefig(['figures/vbmcVofdmVotfsVodss_SymbSNR' simParamStr '.fig']);
            pause(0.1)        
            
            disp(['~~~~~~~~~ Run Time = ' num2str(toc(tExecAlpha)) 's:'...
                ' alphaMax = ' num2str(alphaMax) ','...
                ' tauMax = ' num2str(tauMax*1e3) ' ms,'...
                ' Np = ' num2str(Np) ' paths ~~~~~~~~~']);
        end %end of loop over 'alphaMax'        
                
    end %end of loop over 'tauMax'    
    
end %end of loop over 'Np'
disp(['********* Total run time: ' num2str(toc(tExec)) ' seconds *********']);

%% Plot results on Jain index calculated for VBMC, OFDM and OTFS
figure, 
lgndStr = cell(1,length(tauMaxList)*length(alphaMaxList)); 
iLgnd=0;
for iTau =  [1 length(tauMaxList)] %loop across the tauMax indices to be plotted
    for iAlpha = [1 length(alphaMaxList)]  %loop across the alphaMax indices to be plotted 
        pStr = ['\tau_{max}= ' num2str(tauMaxList(iTau)*1e3) 'ms, \alpha_{max}='...
            num2str(alphaMaxList(iAlpha))];        
        
        plot(NpList, vbmcChanJainAvg(:,iTau,iAlpha),'-o'), hold on,
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['VBMC: ' pStr];        
        plot(NpList, odssChanJainAvg(:,iTau,iAlpha),'-d'), hold on,
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['ODSS: ' pStr];                
        plot(NpList, ofdmChanJainAvg(:,iTau,iAlpha),'-x'),
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['OFDM: ' pStr];        
        plot(NpList, otfsChanJainAvg(:,iTau,iAlpha),'-*'), grid on,
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['OTFS: ' pStr];  
        
        pause(0.01);
    end
end
xlabel('Np'), ylabel('J(h)'), ylim([0 1])
legend(lgndStr{1:iLgnd}), title('Jane''s Index')

%% Plots results of ICI & SIR 
sirIndx = 5;% 1 <= sirInd <= 4+length(betaModn)
for iNp = 1:length(NpList)    
    iNp, %lgndStr = cell(1,length(alphaMaxList)); iLgnd=0;
    for iTau = 1:length(tauMaxList)        
        pStr = ['P=' num2str(NpList(iNp)) ',' ...
            ' \tau_{max}=' num2str(tauMaxList(iTau)*1e3) ' ms']; 
        figure, lgndStr = cell(3,1); iLgnd=0;
        ofdmMaxICIdBAvgiNpiTau = ofdmMaxICIdBAvg(iNp,iTau,:);
        subplot(311),
        semilogx(alphaMaxList-1, ofdmMaxICIdBAvgiNpiTau(:), '-*'),         
        hold on, grid on, 
        ofdmSIRdBAvgiNpiTau = squeeze(ofdmSIRstatsdBAvg(iNp,iTau,:,sirIndx));
        subplot(312),
        semilogx(alphaMaxList-1, ofdmSIRdBAvgiNpiTau, '-*'),    
        hold on, grid on
        ofdmJainIndexAvgiNpiTau = ofdmChanJainAvg(iNp,iTau,:);
        subplot(313),
        semilogx(alphaMaxList-1, ofdmJainIndexAvgiNpiTau(:),'-*'),
        hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['OFDM: ' pStr]; 
        
        otfsMaxICIdBAvgiNpiTau = otfsMaxICIdBAvg(iNp,iTau,:); 
        subplot(311),
        semilogx(alphaMaxList-1, otfsMaxICIdBAvgiNpiTau(:), '-+'),         
        hold on, grid on,         
        otfsSIRdBAvgiNpiTau = squeeze(otfsSIRstatsdBAvg(iNp,iTau,:,sirIndx));
        subplot(312),
        semilogx(alphaMaxList-1, otfsSIRdBAvgiNpiTau, '-+'),  
        hold on, grid on
        otfsJainIndexAvgiNpiTau = otfsChanJainAvg(iNp,iTau,:);
        subplot(313),
        semilogx(alphaMaxList-1, otfsJainIndexAvgiNpiTau(:),'-+'),  
        hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['OTFS: ' pStr]; 
        
        vbmcMaxICIdBAvgiNpiTau = vbmcMaxICIdBAvg(iNp,iTau,:); 
        subplot(311),
        semilogx(alphaMaxList-1, vbmcMaxICIdBAvgiNpiTau(:), '-o'),         
        hold on, grid on,         
        vbmcSIRdBAvgiNpiTau = squeeze(vbmcSIRstatsdBAvg(iNp,iTau,:,sirIndx));
        subplot(312),
        semilogx(alphaMaxList-1, vbmcSIRdBAvgiNpiTau, '-o'),   
        hold on, grid on
        vbmcJainIndexAvgiNpiTau = vbmcChanJainAvg(iNp,iTau,:);
        subplot(313),
        semilogx(alphaMaxList-1, vbmcJainIndexAvgiNpiTau(:),'-o'),   
        hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['VBMC: ' pStr]; 
        
        odssMaxICIdBAvgiNpiTau = odssMaxICIdBAvg(iNp,iTau,:); 
        subplot(311),
        semilogx(alphaMaxList-1, odssMaxICIdBAvgiNpiTau(:), '-o'),         
        hold on, grid on,         
        odssSIRdBAvgiNpiTau = squeeze(odssSIRstatsdBAvg(iNp,iTau,:,sirIndx));
        subplot(312),
        semilogx(alphaMaxList-1, odssSIRdBAvgiNpiTau, '-o'),   
        hold on, grid on
        odssJainIndexAvgiNpiTau = odssChanJainAvg(iNp,iTau,:);
        subplot(313),
        semilogx(alphaMaxList-1, odssJainIndexAvgiNpiTau(:),'-o'),   
        hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['ODSS: ' pStr];         
        
        subplot(311), 
        xlabel('\delta\alpha_{max}'), ylabel('ICI (dB)'), 
        xlim([0.0001, 0.01]), ylim([-15 15]), 
        title('ICI vs \delta\alpha_{max} = \alpha_{max}-1'), 
        %legend(lgndStr{1:iLgnd})
        subplot(312), 
        xlabel('\delta\alpha_{max}'), ylabel('SIR (dB)'), 
        xlim([0.0001, 0.01]), ylim([-15 15]), 
        title('SIR vs \delta\alpha_{max} = \alpha_{max}-1'), 
        legend(lgndStr{1:iLgnd})
        subplot(313), 
        xlabel('\delta\alpha_{max}'), ylabel('J'), 
        xlim([0.0001, 0.01]), %ylim([-15 15]), 
        title('Jain Index vs \delta\alpha_{max} = \alpha_{max}-1'),         
    end
end

%% Plot of SIR
sirIndx = 10;% 1 <= sirInd <= 4+length(betaModn)
for iNp = 1:length(NpList)    
    iNp, %lgndStr = cell(1,length(alphaMaxList)); iLgnd=0;
    for iTau = 1:length(tauMaxList)        
        pStr = ['P=' num2str(NpList(iNp)) ',' ...
            ' \tau_{max}=' num2str(tauMaxList(iTau)*1e3) ' ms']; 
        figure, lgndStr = cell(3,1); iLgnd=0;
%         ofdmMaxICIdBAvgiNpiTau = ofdmMaxICIdBAvg(iNp,iTau,:);
%         subplot(311),
%         semilogx(alphaMaxList-1, ofdmMaxICIdBAvgiNpiTau(:), '-*'),         
%         hold on, grid on, 
        ofdmSIRdBAvgiNpiTau = squeeze(ofdmSIRstatsdBAvg(iNp,iTau,:,sirIndx));
%         subplot(312),
        semilogx(alphaMaxList-1, ofdmSIRdBAvgiNpiTau, '-*'),    
        hold on, grid on
%         ofdmJainIndexAvgiNpiTau = ofdmChanJainAvg(iNp,iTau,:);
%         subplot(313),
%         semilogx(alphaMaxList-1, ofdmJainIndexAvgiNpiTau(:),'-*'),
%         hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['OFDM: ' pStr]; 
        
%         otfsMaxICIdBAvgiNpiTau = otfsMaxICIdBAvg(iNp,iTau,:); 
%         subplot(311),
%         semilogx(alphaMaxList-1, otfsMaxICIdBAvgiNpiTau(:), '-+'),         
%         hold on, grid on,         
        otfsSIRdBAvgiNpiTau = squeeze(otfsSIRstatsdBAvg(iNp,iTau,:,sirIndx));
%         subplot(312),
        semilogx(alphaMaxList-1, otfsSIRdBAvgiNpiTau, '-+'),  
        hold on, grid on
%         otfsJainIndexAvgiNpiTau = otfsChanJainAvg(iNp,iTau,:);
%         subplot(313),
%         semilogx(alphaMaxList-1, otfsJainIndexAvgiNpiTau(:),'-+'),  
%         hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['OTFS: ' pStr]; 
        
%         vbmcMaxICIdBAvgiNpiTau = vbmcMaxICIdBAvg(iNp,iTau,:); 
%         subplot(311),
%         semilogx(alphaMaxList-1, vbmcMaxICIdBAvgiNpiTau(:), '-o'),         
%         hold on, grid on,         
        vbmcSIRdBAvgiNpiTau = squeeze(vbmcSIRstatsdBAvg(iNp,iTau,:,sirIndx));
%         subplot(312),
        semilogx(alphaMaxList-1, vbmcSIRdBAvgiNpiTau, '-o'),   
        hold on, grid on
%         vbmcJainIndexAvgiNpiTau = vbmcChanJainAvg(iNp,iTau,:);
%         subplot(313),
%         semilogx(alphaMaxList-1, vbmcJainIndexAvgiNpiTau(:),'-o'),   
%         hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['VBMC: ' pStr]; 
        
%         odssMaxICIdBAvgiNpiTau = odssMaxICIdBAvg(iNp,iTau,:); 
%         subplot(311),
%         semilogx(alphaMaxList-1, odssMaxICIdBAvgiNpiTau(:), '-o'),         
%         hold on, grid on,         
        odssSIRdBAvgiNpiTau = squeeze(odssSIRstatsdBAvg(iNp,iTau,:,sirIndx));
%         subplot(312),
        semilogx(alphaMaxList-1, odssSIRdBAvgiNpiTau, '-o'),   
        hold on, grid on
%         odssJainIndexAvgiNpiTau = odssChanJainAvg(iNp,iTau,:);
%         subplot(313),
%         semilogx(alphaMaxList-1, odssJainIndexAvgiNpiTau(:),'-o'),   
%         hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['ODSS: ' pStr];         
        
%         subplot(311), 
%         xlabel('\delta\alpha_{max}'), ylabel('ICI (dB)'), 
%         xlim([0.0001, 0.01]), ylim([-15 15]), 
%         title('ICI vs \delta\alpha_{max} = \alpha_{max}-1'), 
        %legend(lgndStr{1:iLgnd})
%         subplot(312), 
        xlabel('\delta\alpha_{max}'), ylabel('SIR (dB)'), 
        xlim([0.0001, 0.01]), %ylim([-10 5]), 
        title('SIR vs \delta\alpha_{max} = \alpha_{max}-1'), 
        legend(lgndStr{1:iLgnd})
%         subplot(313), 
%         xlabel('\delta\alpha_{max}'), ylabel('J'), 
%         xlim([0.0001, 0.01]), %ylim([-15 15]), 
%         title('Jain Index vs \delta\alpha_{max} = \alpha_{max}-1'),         
    end
end

%% Plots results of SINR
sirIndx = 10;% 1 <= sirInd <= 4+length(betaModn)
for iNp = 1:length(NpList)    
    iNp, %lgndStr = cell(1,length(alphaMaxList)); iLgnd=0;
    for iTau = 1:length(tauMaxList)        
        pStr = ['P=' num2str(NpList(iNp)) ',' ...
            ' \tau_{max}=' num2str(tauMaxList(iTau)*1e3) ' ms']; 
        figure, lgndStr = cell(3,1); iLgnd=0;
%         ofdmMaxICIdBAvgiNpiTau = ofdmMaxICIdBAvg(iNp,iTau,:);
%         subplot(311),
%         semilogx(alphaMaxList-1, ofdmMaxICIdBAvgiNpiTau(:), '-*'),         
%         hold on, grid on, 
        ofdmSINRdBAvgiNpiTau = squeeze(ofdmSINRstatsdBAvg(iNp,iTau,:,sirIndx));
%         subplot(312),
        semilogx(alphaMaxList-1, ofdmSINRdBAvgiNpiTau, '-*'),    
        hold on, grid on
%         ofdmJainIndexAvgiNpiTau = ofdmChanJainAvg(iNp,iTau,:);
%         subplot(313),
%         semilogx(alphaMaxList-1, ofdmJainIndexAvgiNpiTau(:),'-*'),
%         hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['OFDM: ' pStr]; 
        
%         otfsMaxICIdBAvgiNpiTau = otfsMaxICIdBAvg(iNp,iTau,:); 
%         subplot(311),
%         semilogx(alphaMaxList-1, otfsMaxICIdBAvgiNpiTau(:), '-+'),         
%         hold on, grid on,         
        otfsSINRdBAvgiNpiTau = squeeze(otfsSINRstatsdBAvg(iNp,iTau,:,sirIndx));
%         subplot(312),
        semilogx(alphaMaxList-1, otfsSINRdBAvgiNpiTau, '-+'),  
        hold on, grid on
%         otfsJainIndexAvgiNpiTau = otfsChanJainAvg(iNp,iTau,:);
%         subplot(313),
%         semilogx(alphaMaxList-1, otfsJainIndexAvgiNpiTau(:),'-+'),  
%         hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['OTFS: ' pStr]; 
        
%         vbmcMaxICIdBAvgiNpiTau = vbmcMaxICIdBAvg(iNp,iTau,:); 
%         subplot(311),
%         semilogx(alphaMaxList-1, vbmcMaxICIdBAvgiNpiTau(:), '-o'),         
%         hold on, grid on,         
        vbmcSINRdBAvgiNpiTau = squeeze(vbmcSINRstatsdBAvg(iNp,iTau,:,sirIndx));
%         subplot(312),
        semilogx(alphaMaxList-1, vbmcSINRdBAvgiNpiTau, '-o'),   
        hold on, grid on
%         vbmcJainIndexAvgiNpiTau = vbmcChanJainAvg(iNp,iTau,:);
%         subplot(313),
%         semilogx(alphaMaxList-1, vbmcJainIndexAvgiNpiTau(:),'-o'),   
%         hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['VBMC: ' pStr]; 
        
%         odssMaxICIdBAvgiNpiTau = odssMaxICIdBAvg(iNp,iTau,:); 
%         subplot(311),
%         semilogx(alphaMaxList-1, odssMaxICIdBAvgiNpiTau(:), '-o'),         
%         hold on, grid on,         
        odssSINRdBAvgiNpiTau = squeeze(odssSINRstatsdBAvg(iNp,iTau,:,sirIndx));
%         subplot(312),
        semilogx(alphaMaxList-1, odssSINRdBAvgiNpiTau, '-o'),   
        hold on, grid on
%         odssJainIndexAvgiNpiTau = odssChanJainAvg(iNp,iTau,:);
%         subplot(313),
%         semilogx(alphaMaxList-1, odssJainIndexAvgiNpiTau(:),'-o'),   
%         hold on, grid on
        iLgnd = iLgnd+1; lgndStr{iLgnd} = ['ODSS: ' pStr]; 
        
%         subplot(311), 
%         xlabel('\delta\alpha_{max}'), ylabel('ICI (dB)'), 
%         xlim([0.0001, 0.01]), ylim([-15 15]), 
%         title('ICI vs \delta\alpha_{max} = \alpha_{max}-1'), 
        %legend(lgndStr{1:iLgnd})
%         subplot(312), 
        xlabel('\delta\alpha_{max}'), ylabel('SINR (dB)'), 
        xlim([0.0001, 0.01]), %ylim([-10 5]), 
        title('SINR vs \delta\alpha_{max} = \alpha_{max}-1'), 
        legend(lgndStr{1:iLgnd})
%         subplot(313), 
%         xlabel('\delta\alpha_{max}'), ylabel('J'), 
%         xlim([0.0001, 0.01]), %ylim([-15 15]), 
%         title('Jain Index vs \delta\alpha_{max} = \alpha_{max}-1'),         
    end
end