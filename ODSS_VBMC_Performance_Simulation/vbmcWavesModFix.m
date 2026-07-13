function [ G, fLn, fHn ] = vbmcWavesModFix( fL, fH, dF, T, Fs, alphaMax, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% B = fH - fL;
betaB = dF*T;
Nsub = floor( log(fH/fL)/log(1 + dF/fL) ) + 1;
f = fL*(1+dF/fL).^(0:Nsub);
% Nsub = floor(log(B/dF*(alphaMax-1)+1)/log(alphaMax));
% f = [0 cumsum(dF*alphaMax.^(0:Nsub-1))]+fL;

fLn = f(1:end-1)*alphaMax; %shrinking the band a bit from the lower edge
fHn = f(2:end)/alphaMax;   %shrinking the band a bit from the upper edge
% figure(1),
% subplot(211), plot(fLn, 'bo'), hold on, grid on,
% plot(fHn, 'r*'), hold off
% subplot(212), plot(fHn-fLn, '>'), grid on
Tmax = T;%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ betaB/min(fHn-fLn);
M = ceil(Fs*Tmax);
G = zeros(M, Nsub);
dF1 = fHn(1) - fLn(1);%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ added during fix
for nSub = 1:Nsub
    if isempty(varargin)
        dFn = fHn(nSub) - fLn(nSub);
        Tn = T*dF1/dFn;%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ betaB/dFn;
        kn = dFn/Tn;
        Mn = ceil(Fs*Tn);
        tn = (0:Mn-1)'/Fs;   
        m1 = floor((M-Mn)/2)+1;
        G(m1:m1+Mn-1,nSub) = ...
            exp(1i*2*pi*(fLn(nSub)*tn + 1/2*kn*tn.^2));
    else
        cwin = varargin{1};
        ct = superChirp(fHn(nSub), fLn(nSub), T, Fs, cwin);
        G(1:length(ct),nSub) = ct;        
    end
end
end