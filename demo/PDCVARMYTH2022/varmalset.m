%% VARMALSET
%        Fit a transfer matrix to multivariate inputs  via Least Squares
%         (limited to identical number of input channels as output channels)
%
%% Syntax
%        [AA,BB,we,pe] = VARMALSET(y,x,p,q)
%
%% Input arguments
%        x      - input 
%        y      - output
%        p      - model order AR part
%        q      - model order (q+1 output matrices)
%
%% Output arguments
%        AA     - [m,m,p] array of AR parameters 
%        BB     - [m,m,q+1] array of parameters
%        we     - model observation errors
%        pe     - model observation error covariance matrix
%

function [AA,BB,we,pe]=varmalset(y,x,p,q)

[m,n] = size(y);
yo = y';

xi = x';

HB = [];HA = HB;
for i = 1:m
   VA = toeplitz(yo(:,i),[yo(1,i) zeros(1,p)]);
   VA = VA(:,2:end);
   HA = [HA VA];
end

GA = zeros(m*n,m*p);
for i = 1:m
   GA(1+(i-1)*n:i*n, 1+(p)*m*(i-1):i*(p)*m) = HA;
end

for i = 1:m
   VB = toeplitz(xi(:,i),[xi(1,i) zeros(1,q)]);
   VB = VB(:,2:end);
   HB = [HB VB];
end

GB = zeros(m*n,m*(q));
for i = 1:m
   GB(1+(i-1)*n:i*n, 1+(q)*m*(i-1):i*(q)*m) = HB;
end

G = [GA GB];
T = G'*(yo(:)-xi(:));
B = inv(G'*G) * T;
we = yo(:)-G*B; we = reshape(we,n,m)'; pe = we*we'/n;
v = B(1:m*m*p);
vv = reshape(v,p,m*m);

k = 1;
for i = 1:m,
   for j = 1:m,
      AA(i,j,:) = vv(:,k);
      k = k+1;
   end
end

vt = B(m*m*p+1:end);
vvt = reshape(vt,q,m*m);
k = 1;
for i = 1:m
   for j = 1:m
      BB(i,j,:) = vvt(:,k);
      k = k+1;
   end
end
