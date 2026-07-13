function [SINRdB, sinrJain, gamma] = ... %[minSIRdB, maxSIRdB, meanSIRdB, medianSIRdB, varargout] =...
    SINRstats(H, noisePower, varargin)
%ICI.m Summary of this function goes here
%   Detailed explanation goes here

gamma = (abs(diag(H)).^2)./(sum(abs( H - diag(diag(H)) ).^2, 2)+ noisePower);%gamma: Nx1 vector, H:NxN matrix

SINRdB = real([min(10*log10( gamma )) max(10*log10( gamma ))...
    10*log10( mean(gamma) ) 10*log10( median(gamma) )]);
if ~isempty(varargin)
    beta = varargin{1};
    eeSIR = beta(:)'.*log(mean(exp(-gamma(:)*(1./beta(:)'))));
    SINRdB = [SINRdB real(10*log10(eeSIR))];
end

avgGamma = sum(abs(diag(H)).^2)/sum(sum(abs( H - diag(diag(H)) ).^2, 2)+...
    noisePower);
SINRdB = [SINRdB avgGamma];

sinrJain = JainIndex(gamma);
end