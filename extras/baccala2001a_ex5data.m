function u=baccala2001example5data(nPoints)
% Baccala & Sameshima. Partial directed coherence: a new concept in neural
% structure determination. Biol. Cybern. 84:463-474, 2001.
%                http://dx.doi.org/10.1007/PL00007990
% 
% Example VAR(2) five-variable with loop and feedback

N=10000;
disp('======================================================================');
disp('             Linear VAR(2) Model Example 5 - closed-loop')
disp('         Baccala & Sameshima. Biol. Cybern. 84:463-474, 2001.')
disp('           x1==>x2  x2-->x3 x3-->x4 x4-->x5 x5-->x4 x5-->x1');
disp('======================================================================');

randn('state', sum(100*clock))
wi=randn(5,N);
x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);
x4=zeros(1,N);
x5=zeros(1,N);
% Variables initialization
for t=1:4,
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1); x4(t)=randn(1);
   x5(t)=randn(1);
end;

for t=5:N,
   x1(t) = 0.95*sqrt(2)*x1(t-1) - 0.9025*x1(t-2) + 0.5*x5(t-2)+ wi(1,t);
   x2(t) = -0.5*x1(t-1) + wi(2,t);
   x3(t) = 0.4*x2(t-2) + wi(3,t);
   x4(t) = -0.5*x3(t-1) + 0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) + wi(4,t);
   x5(t) = -0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) + wi(5,t);
end;

nDiscard=1000; % number of points discarded at beginning of series
%nPoints       % number of analyzed samples points
if (nDiscard+nPoints) > N, error('Too many data points are been requested. Edit baccala2001example5data.m file.'); end;

y=[x1' x2' x3' x4' x5']; % data must be organized column-wise
u=y(nDiscard+1:nDiscard+nPoints,:);
