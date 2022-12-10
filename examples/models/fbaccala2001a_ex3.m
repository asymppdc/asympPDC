%% FBACCALA2001A_EX3
% Linear five-dimensional VAR[3] with feedback model from Baccala & Sameshima.
% Biol.Cybern.84:463-74, 2001. Example 3 (pag.468)') with following interactions
% x1-->x2  x1-->x3 x1-->x4 x4-->x5 x5-->x4');
%                http://dx.doi.org/10.1007/PL00007990

% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.


function [ u ] = fbaccala2001a_ex3( nPoints, nBurnIn )

disp('======================================================================');
disp('               Linear five-dimensional VAR[3] Model')
disp(' Baccala & Sameshima. Biol.Cybern.84:463-74, 2001. Example 3 (pag.468)')
disp('             x1-->x2  x1-->x3 x1-->x4 x4-->x5 x5-->x4');
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

N=nBurnIn+nPoints; % Number of data points to be generated.

% Variables initialization.

rng('default')
rng('shuffle')

ei=randn(5,N);
x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);
x4=zeros(1,N);
x5=zeros(1,N);

for t=1:4
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1); x4(t)=randn(1);
   x5(t)=randn(1);
end

for t=5:N
   x1(t) = 0.95*sqrt(2)*x1(t-1) - 0.9025*x1(t-2) + ei(1,t);
   x2(t) = 0.5*x1(t-2) + ei(2,t);
   x3(t) = -0.4*x1(t-3) + ei(3,t);
   x4(t) = -0.5*x1(t-2) + 0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) + ei(4,t);
   x5(t) = -0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) + ei(5,t);
end

y=[x1' x2' x3' x4' x5']; % data must be organized column-wise

u=y(nBurnIn+1:N,:);

end

