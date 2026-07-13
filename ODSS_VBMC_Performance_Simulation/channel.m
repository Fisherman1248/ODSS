function [ r, Ps, Pn, rs, w ] = channel( s, fs, taup, alphap, hp, SNR, varargin )
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Authors : (1) Arunkumar K. P.
%          (2) Chandra R. Murthy
%Address : (1) Ph.D. Scholar,
%              Signal Processing for Communications Lab, ECE Department,
%              Indian Institute of Science, Bangalore, India-560 012.
%          (2) Professor,
%              Electrical Communication Engineering (ECE) Department,
%              Indian Institute of Science, Bangalore, India-560 012.
%Email   : arunkumar@iisc.ac.in
%
%Revision History
% Version : 2.0
% Revisions: 12-11-2020
% Last Revision: 28-02-2021
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% This is a script under development. The work on ODSS communication scheme
% is planned for a joint patent filing by the authors through institute
% (IISc)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%   Detailed explanation goes here
        if length(nargin)>1
            fsR = varargin{2};
            s = resample(s, fsR/fs, 1);
        else
            fsR = fs;
        end
        Np = length(taup); %number of paths
        % (1) Obtain signal after passing s(t) through the simulated delay-scale channel
        rs = [];%signal after passing s(t) through the simulated delay-scale channel
        for nPath = 1:Np %add time-scaled and delayed arrivals via each path
            rn = [zeros(1,round(taup(nPath)*fsR))...
                resample(s,round(alphap(nPath)*1e4),10000)...
                ]; %scaled and delayed version of the transmitted signal s(t)
            lenMax = max(length(rs),length(rn));%to make rn and rs to same length vectors
            rs(end+1:lenMax) = 0;%make rn and rs to same length vectors
            rn(end+1:lenMax) = 0;%make rn and rs to same length vectors
            rs = rs + hp(nPath)*rn;%add scaled and delayed version to the received signal
        end
        
        if length(nargin)>1
            fsR = varargin{2};
            rs = resample(rs, 1, fsR/fs);
        end        
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % (2) Add receiver noise (AWGN)
        Ps = mean( abs(rs).^2 );%received signal (waveform level) power
        Ps0 = Np; %true or expected signal power
        Pn = Ps0*10^(-SNR/10); %noise (waveform level) power for the specified SNR
        w = sqrt(Pn)*(randn(size(rs))+1i*randn(size(rs)))/sqrt(2); %AWGN
        r = rs + w; %received waveform after (additive) noise corruption        
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        if ~isempty(varargin)
            trimStr = varargin{1};
            if strcmpi(trimStr, 'trim')
                if length(r) >= length(s)
                    r = r(1:length(s));
                else
                    r = [r zeros(1,length(s) - length(r))];
                end
            end
        end
%         r = r(1:length(s));
        
end
