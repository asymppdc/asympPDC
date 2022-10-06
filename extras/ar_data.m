function c = ar_data(A, er, m, ndiscard)
% function c = ar_data(A, er, m, ndiscard)
%    '''Simulate ar-model from A matrix
%
%       Input:
%         A(n, n, p) - AR model (n - number of signals, p - model order)
%         er(n) - variance of innovations
%         m - length of simulated time-series
%
%       Output:
%         data(n, m) - simulated time-series
%     '''

% (C) Koichi Sameshima & Luiz A. BaccalÁ, 2022. 
% See file license.txt in installation directory for licensing terms.

if ndims(A) == 2
   [n,~]=size(A);
   p=1;
elseif ndims(A) == 3
   [n,~,p] = size(A);
else
   error('A matrix dimension > 3.');
end

if isempty(er)
   er = eye(n); %identity(n)
end
if min(size(er))==1
   er = diag(er);
end
randn('state', sum(100*clock));
er = diag(er);

w = randn(n,m+ndiscard+p);

for kk=1:m+ndiscard-p, w(:,kk)=er.*w(:,kk); end

data = zeros(n, m+ndiscard);
for i =p+1:m+ndiscard
   for j=1:p
       data(:,i) = data(:,i) + A(:,:,j)*data(:,i-j);
   end
   data(:,i) = data(:,i) + w(:,i);
end

c = data(:,ndiscard+1:m+ndiscard);
