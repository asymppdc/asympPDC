function [ u ] = fbaccala2001b_model1_feedback( nPoints, nDiscard )

% Baccala & Sameshima. Overcoming the limitations of correlation analysis 
% for many simultaneously processed neural structures, Progress in Brain 
% Research, 130:33--47, 2001.
%            http://dx.doi.org/10.1016/S0079-6123(01)30004-3
% 
% Example Model I - 7-dimensional VAR[2] model with loop and feedback


disp('======================================================================');
disp('                       Linear VAR[2] Model I')
disp('        Baccala & Sameshima. Prog Brain Research, 130:33--47, 2001.')
disp('  x1-->x2  x1-->x3  x2-->x3 x3-->x4 x4==>x5 x5-->x1  x5-->x4 x6==>x7');
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
ei=randn(7,N);
x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);
x4=zeros(1,N);
x5=zeros(1,N);
x6=zeros(1,N);
x7=zeros(1,N);
% Variables initialization
for t=1:4,
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1); x4(t)=randn(1);
   x5(t)=randn(1); x6(t)=randn(1); x7(t)=randn(1); 
end;

for t=5:N,
   x1(t) = 0.95*sqrt(2)*x1(t-1) - 0.9025*x1(t-2) + 0.5*x5(t-2) + ei(1,t);
   x2(t) = -0.5*x1(t-1) + ei(2,t);
   x3(t) =  0.2*x1(t-1) + 0.4*x2(t-2) + ei(3,t);
   x4(t) = -0.5*x3(t-1) + 0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-2) + ei(4,t);
   x5(t) = -0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) + ei(5,t);
   x6(t) = 0.95*sqrt(2)*x6(t-1) - 0.9025*x6(t-2) + ei(6,t);
   x7(t) = -0.1*x6(t-2) + ei(7,t);   
end;

y=[x1' x2' x3' x4' x5' x6' x7']; % data must be organized column-wise
u=y(nDiscard+1:N,:);

end

