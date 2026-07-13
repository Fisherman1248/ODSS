function [y] = delay(x,tau,fs,tauMax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
L = ceil(tauMax*fs);
xz = [x zeros(1, L)];
M = length(xz);
Xz = fft(xz);
Y = Xz.*exp(-1i*2*pi*(0:M-1)*fs/M*tau);
y = ifft(Y, 'symmetric');
% y = y(1:length(x));
end

