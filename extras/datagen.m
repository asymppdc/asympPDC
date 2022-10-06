%% DATAGEN
%        Generate multivariate time series from autoregressive coefficients A  and
%        covariance matrix, pf.
% 
%% Syntax:
%         [y,seed_out]=DATAGEN(A,pf,Ndata,R,seed)
%
%% Input Arguments: 
%   A:          NxNxIP system definition matrix, IP is model order
%   pf:         NxN Input covariance matrix
%   Ndata:      Number of time samples
%   R:          Number of realizations (default = 1)
%   seed:       Seed default - internally defined (nargin<5)
%        
%   nBurnIn=10000 (internally defined)
%
%% Output Argument:
%   y - N x Time (=Ndata) x R

% July 09, 2021, LAB
%

% future seed inclusion
% default burn-in period to deal with transients

% (C) Koichi Sameshima & Luiz A. BaccalÁ, 2022. 
% See file license.txt in installation directory for licensing terms.

%%


function [y,seed_out]=datagen(A,pf,Ndata,R,seed)

nBurnIn=10000;

%% Warning on rng usage
%    Error using rng (line xx)
% The current random number generator is the legacy generator.  This is
% because you have executed a command such as rand('state',0), which
% activates MATLAB's legacy random number behavior.  You may not use RNG to
% reseed the legacy random number generator.
% 
% Use rng('default') to reinitialize the random number generator to its
% startup configuration, or call RNG using a specific generator type, such
% as rng(seed,'twister').
if nargin==5
   rng('default')
   rng(seed)

end
seed_out=rng;

if nargin<4
   R=1;
end
[n,~,p]=size(A);
V=chol(pf);
epsilon=randn(n,Ndata + nBurnIn + p,R);
y=zeros(size(epsilon));
for r=1:R
   epsilon(:,:,r)=V'*epsilon(:,:,r);
end
for r=1:R
   u=zeros(n,Ndata + nBurnIn + p);v=epsilon(:,:,r);
   for t=p + 1:Ndata + nBurnIn + p
      ut=zeros(n,1);
      for i=1:p
         ut=A(:,:,i)*u(:,t-i) + ut;
      end
      u(:,t)=ut;
      u(:,t)=u(:,t) + v(:,t);
   end
   y(:,:,r)=u(:,:);
end
y=y(:,nBurnIn + p + 1:Ndata + nBurnIn + p,:);
end

