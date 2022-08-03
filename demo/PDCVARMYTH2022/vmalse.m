%% VMALSE
%       Fit a transfer matrix to multivariate inputs via Least Squares
%       (limited to identical number of input channels as output channels)
%
%% Syntax
%       [BB,we,pe]=VMALSE(y,x,q)
%
%% Input arguments
%       x       - input 
%       y       - output
%       q       - model order (q+1 output matrices)
%
%% Output arguments 
%       BB      - [m,m,q+1] array of parameters
%       we      - model observation errors
%       pe      - model observation error covariance matrix
%

%       B       - auxiliary  BB output -early version)

function [BB,we,pe]=vmalse(y,x,q);

[m,n] = size(y);
yo = y'; yo = yo(:);
xi = x';

H = [];
for i = 1:m
    H = [H toeplitz(xi(:,i),[xi(1,i) zeros(1,q)])];
end

G = zeros(m*n,m*(q+1));
for i = 1:m
    G(1+(i-1)*n:i*n,1+(q+1)*m*(i-1):i*(q+1)*m) = H;
end
T = G'*yo;
B = inv(G'*G)*T;

we = yo-G*B; we = reshape(we,n,m)'; pe = we*we'/n;

B = reshape(B,q+1,m*m)';
for i = 1:q+1
    BB(:,:,i) = reshape(B(:,i),m,m)';
end


