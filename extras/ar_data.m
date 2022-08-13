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

if ndims(A) == 2,
   [n,~]=size(A);
   p=1;
elseif ndims(A) == 3,
   [n,~,p] = size(A);
else
   error('A matrix dimension > 3.');
end;

if isempty(er),
   er = eye(n); %identity(n)
end;
if min(size(er))==1,
   er = diag(er);
end;
randn('state', sum(100*clock));
er = diag(er);
%w = mnorm(zeros(n), er, m+ndiscard-p)
w = randn(n,m+ndiscard+p);
%w = er.*randn(n,n,m+ndiscard-p);

for kk=1:m+ndiscard-p, w(:,kk)=er.*w(:,kk); end;

data = zeros(n, m+ndiscard);
for i =p+1:m+ndiscard,
   for j=1:p,
%      data(:,i) = data(:,i) + dot(A(:,:,j), data(:,i-j-1)); % checar os indices
       data(:,i) = data(:,i) + A(:,:,j)*data(:,i-j); % checar os indices
   end;
   data(:,i) = data(:,i) + w(:,i);
end;
%print time.clock()

c = data(:,ndiscard+1:m+ndiscard);

%save

% function A = ar_models(id, lam),
% %function A = ar_models(id, lam = 0.0),
% 
% if nargin < 2, lam = 0.; end;
% 
% models = [
%    %0
%    [array([[[0.2, 0],[0, 0],[0.3,-0.2]],
%    [[0, 0],[0.8,-0.1],[0.4,-0.1]],
%    [[0, 0],[-0.1,0.2],[0.4,0.1]]], dtype = float),
%    identity(3)],
%    %1
%    [array([[[4,-4],[3,3]],[[0,-2],[2,-3]]], dtype=float).reshape(2,2,2)/20,
%    array([[0.7,0],[0,2]], dtype = float)],
%    %2 sunspot melanoma
%    sun,
%    %3
%    [array([[[4,3,-2],[-2,-5,3]],[[4,-2,1],[-4,0,3]]], dtype=float).reshape(2,2,3)/20,
%    array([[0.7,0.3],[0.3,2]], dtype = float)],
%    %4 JAS Daniel (12)
%    [array([[[0.2, 0],[-0.4, -0.2],[0.3,0]],
%    [[lam, 0],[0.8,-0.1],[0.4,0.0]],
%    [[0, 0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float),
%    identity(3)],
%    %5 Stein 2010
%    [array([[[0.2, 0],[-0.04, -0.02],[0.3,0]],
%    [[lam, 0],[0.8,-0.1],[0.4,0.0]],
%    [[0, 0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float),
%    diag(array([0.1, 5, 2], dtype = float))],
% 
%    ]
% A = models[id];