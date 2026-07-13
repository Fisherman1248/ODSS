function [maxICIdB] = maxICI(H)
%ICI.m Summary of this function goes here
%   Detailed explanation goes here

maxICIdB =... 
    max(10*log10(sum(abs(H - diag(H)),2))./max(abs(diag(H))));
%     max( max( 20*log10(bsxfun(@times, abs(H - diag(diag(H))),...
%     1./diag(abs(H))')) ) );
end
