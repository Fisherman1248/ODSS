function [JainX] = JainIndex(x)
%JainIndex.m Summary of this function goes here
%   Computes the Raj Jain Index of the vector x defined by:
% J(x) = abs(x.'*1n).^2/(n*x'*x), where n = length(x)
%
n = length(x);
JainX = abs(sum(x))^2/(n*(x(:)'*x(:)));
end
