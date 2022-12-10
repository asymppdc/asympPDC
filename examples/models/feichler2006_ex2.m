function [ u ] = feichler2006_ex2( nPoints, nBurnIn, flgManual)


% Eichler. On the evaluation of information flow in multivariate systems
% by the directed transfer function.
%            Biol Cybern (2006) 94: 469?482
%
% <http://dx.doi.org/10.1007/s00422-006-0062-z>
%
%  Example - two-dimensional VAR[4].

% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.

disp('======================================================================');
disp('                   Two-dimensional linear VAR[4] Model')
disp('                           Eichler (2006).')
disp('                               x1-->x2');
disp('======================================================================');

if (nargin == 0)
   nPoints = 1000;
   nBurnIn = 5000;
   disp(['Adopting default ' int2str(nBurnIn) ' discarding points, and'])
   disp(['generating ' int2str(nPoints) ' simulation data point.'])
   
elseif   (nargin < 2)
   nBurnIn = 5000;
   disp(['Adopting default ' int2str(nBurnIn) ' discarding points.'])
end

if (nBurnIn < 1)
   nBurnIn = 5000;
   disp(['Adopting default ' int2str(nBurnIn) ' discarding points.'])
end

if nPoints < 10
   nPoints = 100;
   disp(['Adopting default ' int2str(nPoints) ' simulation data points.'])
end

N = nBurnIn+nPoints; % number of simulated points.


if ~exist('flgManual','var')
   flgManual = 0;
end

if flgManual
   if exist('randn_manual_state.mat','file')
      load randn_manual_state
      randn('state', s);
      disp(['Using state saved in "randn_manual_state.mat" to reproduce figure in the manual.']);
   else
      disp(['Did not find "randn_manual_state.mat" file.' ...
         'Assigned "sum(100*clock)" initial state.']);
      randn('state', sum(100*clock))
   end
else
   rng('shuffle')
   disp('Seeds the random number generator using rng(''shuffle'').')
end
% Variables initialization
e=randn(2,N);
x1=zeros(1,N);
x2=zeros(1,N);

disp('======================================================================');

for t=1:4
   x1(t)=randn(1); x2(t)=randn(1);
end

chLabels = {'X_1';'X_2'}; %or %chLabels = [];

for t=5:N
   x1(t) = 7/8*x1(t-1) - 4/5*x1(t-2) +4/5*x1(t-3) - 1/2*x1(t-4) + e(1,t);
   x2(t) = 1/20*x1(t-1) + 1/20*x1(t-2) + e(2,t);
   
end

y=[x1' x2']; % data must be organized column-wise
u=y(nBurnIn+1:N,:);
