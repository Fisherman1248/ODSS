function [SIRdB] = SIR(H)
%ICI.m Summary of this function goes here
%   Detailed explanation goes here

SIRdB =... 
    min(10*log10( (abs(diag(H)).^2)./sum(abs( H - diag(diag(H)) ).^2, 2) ));
end