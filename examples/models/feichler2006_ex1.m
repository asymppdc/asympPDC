function [u ] = feichler2006_ex1 (nPoints, nBurnIn, flgManual)

% Eichler. On the evaluation of information flow in multivariate systems 
% by the directed transfer function.
%            Biol Cybern (2006) 94: 469?482
%
% <http://dx.doi.org/10.1007/s00422-006-0062-z> 
% 
%  Example - Three-dimensional VAR[2].

% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.


disp('======================================================================');
disp('                   Three dimensional linear VAR[3] Model')
disp('                           Eichler, 2006.')
disp('                      x1-->x2  x2-->x1 x2-->x3 ');
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
e=randn(3,N);
x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);

disp('======================================================================');

for t=1:3
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1);
end

for t=4:N
   x1(t) = 4/15*x1(t-1) - (1./4.)*x1(t-2) + .4*x2(t-1) - .2*x2(t-2) + e(1,t);
   x2(t) = 0.4*x1(t-1) + 0.2*x1(t-2) + (5./28.)*x2(t-1) - (1./9.)*x2(t-2) + e(2,t);
   x3(t) = 0.4*x2(t-1) + 0.2*x2(t-2) + (5./28.)*x3(t-1) - (1./9.)*x3(t-2) + e(3,t);
end

y = [x1' x2' x3']; 
u = y(nBurnIn+1:N,:); % data must be organized row-wise

%end;