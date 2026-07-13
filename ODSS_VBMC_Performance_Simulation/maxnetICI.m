function [maxICIdB] = maxnetICI(H)
%ICI.m Summary of this function goes here
%   Detailed explanation goes here

maxICIdB =... 
    max(10*log10( sum(abs(H - diag(diag(H))).^2,2)./(abs(diag(H)).^2) ));
end