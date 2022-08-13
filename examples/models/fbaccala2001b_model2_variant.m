function [ u ] = fbaccala2001b_model2_variant( nPoints, nDiscard )

% Variation of model II from
% Baccala & Sameshima. Overcoming the limitations of correlation analysis 
% for many simultaneously processed neural structures, Progress in Brain 
% Research, 130:33--47, 2001.
%            http://dx.doi.org/10.1016/S0079-6123(01)30004-3
% 
% Example 6-dimensional VAR[4] Scalp EEG model II *variant*

% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.


disp('======================================================================');
disp('                       Linear VAR[4] Variant Model II')
disp('        Baccala & Sameshima. Prog Brain Research, 130:33--47, 2001.')
disp('        x1==>x2  x1==>x4  x2==>x3 x2-->x5 x5==>x6 x6-->x4  x6==>x5');
disp('======================================================================');

if (nargin == 0),
   nPoints = 1000;
   nDiscard = 5000;
   disp(['Adopting default ' int2str(nDiscard) ' discarding points, and'])
   disp(['generating ' int2str(nPoints) ' simulation data point.'])
   
elseif   (nargin < 2),
   nDiscard = 5000;
   disp(['Adopting default ' int2str(nDiscard) ' discarding points.'])
end;

if (nDiscard < 1),
   nDiscard = 5000;
   disp(['Adopting default ' int2str(nDiscard) ' discarding points.'])
end;

if nPoints < 10,
   nPoints = 100;
   disp(['Adopting default ' int2str(nPoints) ' simulation data points.'])
end;

N = nDiscard+nPoints; % number of simulated points.


randn('state', sum(100*clock))
ei=randn(6,N);
x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);
x4=zeros(1,N);
x5=zeros(1,N);
x6=zeros(1,N);

% Variables initialization
for t=1:4,
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1); x4(t)=randn(1);
   x5(t)=randn(1); x6(t)=randn(1);; 
end;

for t=5:N,
   x1(t) = 0.95*sqrt(2)*x1(t-1) - 0.9025*x1(t-2)+ ei(1,t); % modified w oscillation
%  x1(t) = 1.8982*x1(t-1) - 0.9025*x1(t-2)+ ei(1,t); % original equation
   x2(t) = 0.9*x1(t-2) + ei(2,t);
   x3(t) = 0.85*x2(t-2) + ei(3,t);
   x4(t) = 0.82*x1(t-2) + 0.6*x6(t-3) + ei(4,t);
   x5(t) = -0.9*x6(t-2) + 0.4*x2(t-4) + ei(5,t);
   x6(t) = 0.9*x5(t-2) + ei(6,t);
end;

y=[x1' x2' x3' x4' x5' x6']; % data must be organized column-wise
u=y(nDiscard+1:N,:);

end

