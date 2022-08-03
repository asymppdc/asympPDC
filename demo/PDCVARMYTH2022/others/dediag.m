%% DEDIAG
%       Calculate the correlation matrix from the covariance matrix
%
%% Syntax
%       [Y,x]=DEDIAG(S,invflag)
%
%% Input arguments
%       S       - (N x N) covariance matrix
%       invflag - If any value is provided, calculate partial correlations.
%
%% Output arguments 
%       x       - root square variances
%       Y       - correlation matrix (partial correlations if invflag=1)
%

function [Y,x]=dediag(S,invflag)

x  = sqrt(diag(S));
ix = diag(1./x);
Y  = ix * S * ix;

if nargin > 1
   Y  = inv(Y);
   y  = sqrt(diag(Y));
   iy = diag(1./y); 
   Y  = iy * Y * iy;
   x  = y;
end
