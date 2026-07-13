function [Hsinc] = delayScaleChannelBasis(hp, taup, alphap, W, fc, Fs, Ntx, tauMax)
%delayScaleChannelBasis.m: Summary of this function goes here
% Constructs the time domain baseband channel matrix, Hsinc, for a delay
% and scale spread channel with Np paths each associated with a path delay
% of taup, time-scale alphap, and complex path amplitude hp. The columns of
% Hsinc are the path modified versions of the columns of the Bsinc matrix
% whose columns in turn are the basis functions for baseband signals. We
% have: Bsinc(t,n) = sinc(W(t-nTs)), t=(0:N-1)/Ts. Therefore,
% Hsinc(t,n)= \sum_p hp Bsinc(alphap (t - taup),n)
%----------------------------------------------------------------------
Np = length(hp);%number of paths
% Ntx = ceil(Ttx*Fs); %number samples in the transmitted waveform
Ttx = Ntx/Fs; %transmitted waveform duration
Trx = Ttx + tauMax;  %max. duration of the received signal, in s: Trx = T + tauMax
Nrx = ceil(Trx*Fs/min(alphap))+1; %(3) N = number of channel output time samples observed at the receiver front-end in time-domain
t = (0:Nrx-1)'/Fs;%sampling times
Hsinc = zeros(Nrx,Ntx); %delay-scale channel basis matrix
for p = 1:Np
    Tp = Ttx/alphap(p); %duration of the time-scaled signal arriving along the p'th path
    n1p = ceil(taup(p)*Fs)+1; %sample index relative to first arrival from when the p'th path response starts at the receiver
    n2p = n1p + floor(Tp*Fs);%floor((Tp+taup(p))*Fs)+1; %sample index relative to first arrival when the p'th path response end at the receiver
    tp = alphap(p)*(t(n1p:n2p)-taup(p));
    Hsinc(n1p:n2p,:) = Hsinc(n1p:n2p,:) + ...
        hp(p)*sqrt(alphap(p))*sinc(W* (tp - (0:Ntx-1)/Fs )).*exp(1i*2*pi*fc*tp).*...
        exp(-1i*2*pi*fc*t(n1p:n2p));
end
end