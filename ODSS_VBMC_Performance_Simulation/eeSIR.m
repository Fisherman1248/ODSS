function [eeSIRdB] = eeSIR(H,beta)
%eeSIR.m Summary of this function goes here
%   Detailed explanation goes here
gamma = (abs(diag(H)).^2)./sum(abs( H - diag(diag(H)) ).^2, 2);%gamma: Nx1 vector, H:NxN matrix
eeSIR = beta(:)'.*log(mean(exp(-gamma(:)*(1./beta(:)'))));
eeSIRdB = 10*log10(eeSIR);
end