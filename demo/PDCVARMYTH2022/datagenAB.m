%% DATAGENAB
%       Generate data sample with model parameters provided by A and B matrices.
%
%% Syntax 
%       [y,seed_out,epsilon0]=DATAGENAB(A,B,pf,Ndata,R,seed)
%
%% Input arguments
%       A        - (nChannels X nChannels x IP) system definition matrix
%       B        - (nChannels x nChannels x IQ+1) FIR definition matrices B(:,:,1)=eye(nChannels) - see below
%       pf       - input covariance matrix
%       Ndata    - number of time samples
%       R        - number of realizations (default = 1)
%       seed     - seed default, internally defined if nargin < 6
%
%% Output Arguments: 
%       y        - (nChannels x Ndata x R) data sample 
%       seed_out - used seed
%       epsilon0 - Innovation
%

% 09/07/2021
%

function [y,seed_out,epsilon0]=datagenAB(A,B,pf,Ndata,R,seed)

% future seed inclusion
% Default burnin period to deal with transients
BurnIn = 10000; 

if nargin==6
   % Initialize rng to avoid predictability of sessions
   if isOctave()
      randn('seed',seed); % set randn initial 'seed' randomization.
      seed_out = seed;
   else
      rng(seed);
      seed_out = rng;
   end
end

if nargin < 6
   if isOctave()
      seed_out = sum(100*clock);
      randn('seed',seed_out); % set randn initial 'seed' randomization.
   else
      rng('default')
      rng('shuffle')
      seed_out = rng;
   end
end

if nargin<5
   R = 1;
end

[nChA,~,p] = size(A);
if nChA == 0
   disp('This is a model with A = [].');
end

[nChannels,~] = size(pf);

V = chol(pf);
epsilon = randn(nChannels,Ndata+BurnIn+p,R);
y = zeros(size(epsilon));

for r = 1:R
   epsilon(:,:,r) = V'*epsilon(:,:,r);
end
epsilon0 = epsilon(:,BurnIn+p+1:Ndata+BurnIn+p,:);

if ~isempty(B)
   epsilon = mfir(epsilon,B);
end

if ~isempty(A)
   for r = 1:R
      u = zeros(nChannels,Ndata+BurnIn+p);v = epsilon(:,:,r);
      for t = p+1:Ndata+BurnIn+p
         ut = zeros(nChannels,1);
         for i = 1:p
            ut = A(:,:,i)*u(:,t-i)+ut;
         end
         u(:,t) = ut;
         u(:,t) = u(:,t)+v(:,t);
      end
      y(:,:,r) = u(:,:);
   end
else
   y = epsilon;
end
y = y(:,BurnIn+p+1:Ndata+BurnIn+p,:);
