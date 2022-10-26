function [ u ] = fbaccala2001a_ex4( nPoints, nBurnIn )

% Baccala & Sameshima. Partial directed coherence: a new concept in neural
% structure determination. Biol. Cybern. 84:463-474, 2001.
%                http://dx.doi.org/10.1007/PL00007990
% 
% Five-dimensional VAR[2] without feedback example.

disp('======================================================================');
disp('            Five-dimensional linear VAR[2] Model Example 4')
disp('           Baccala & Sameshima. Biol. Cybern. 84:463-474, 2001.')
disp('               x1==>x2  x2-->x3 x3-->x4 x4-->x5 x5-->x4');
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
% 
% clear; clc
% 
% nBurnIn = 1000;    % number of points discarded at beginning of simulation
% nPoints  = 2000;   % number of analyzed samples points
N=nBurnIn+nPoints; %

randn('state', sum(100*clock))
ei=randn(5,N);
x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);
x4=zeros(1,N);
x5=zeros(1,N);
% Variables initialization
for t=1:4
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1); x4(t)=randn(1);
   x5(t)=randn(1);
end

chLabels = []; % or 
%chLabels = {'x_1';'x_2';'x_3';'x_4';'x_5'};

for t=5:N
   x1(t) = 0.95*sqrt(2)*x1(t-1) - 0.9025*x1(t-2) + ei(1,t);
   x2(t) = -0.5*x1(t-1) + ei(2,t);
   x3(t) = 0.4*x2(t-2) + ei(3,t);
   x4(t) = -0.5*x3(t-1) + 0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) + ei(4,t);
   x5(t) = -0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) + ei(5,t);
end

y=[x1' x2' x3' x4' x5']; % data must be organized row-wise
u=y(nBurnIn+1:N,:)';

end

