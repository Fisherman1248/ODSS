function [maxICIdB] = ICI(H)
%ICI.m Summary of this function goes here
%   Detailed explanation goes here

maxICIdB =...
    max( max( 20*log10(bsxfun(@times, abs(H - diag(diag(H))),...
    1./diag(abs(H))')) ) );
end
