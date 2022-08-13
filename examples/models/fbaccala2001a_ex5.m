function [ u ] = fbaccala2001a_ex4( nPoints, nDiscard )

%
% Baccala & Sameshima. Partial directed coherence: a new concept in neural
% structure determination. Biol. Cybern. 84:463-474, 2001.
%                http://dx.doi.org/10.1007/PL00007990
% 
% Example Five-dimensional VAR[2] with loop and feedback

% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.

disp('======================================================================');
disp('        Pentavariate linear VAR[2] Model Example 5 - closed-loop')
disp('         Baccala & Sameshima. Biol. Cybern. 84:463-474, 2001.')
disp('           x1==>x2  x2-->x3 x3-->x4 x4-->x5 x5-->x4 x5-->x1');
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

%
randn('state', sum(100*clock))

%randn('state', 100*pi); % This command line was used to repeat 
%                        % the analysis with the same data set to check 
%                        % alg and criterion selection.

%% Variables initialization
wi=randn(5,N); x1=zeros(1,N); x2=x1; x3=x1; x4=x1; x5=x1;

for t=1:4,
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1); x4(t)=randn(1);
   x5(t)=randn(1);
end;

chLabels = []; % or  = {'x_1';'x_2';'x_3';'x_4';'x_5'};

for t=5:N,
   x1(t) = 0.95*sqrt(2)*x1(t-1) - 0.9025*x1(t-2) + 0.5*x5(t-2)+ wi(1,t);
   x2(t) = -0.5*x1(t-1) + wi(2,t);
   x3(t) = 0.4*x2(t-2) + wi(3,t);
   x4(t) = -0.5*x3(t-1) + 0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) + wi(4,t);
   x5(t) = -0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) + wi(5,t);
end;

y=[x1' x2' x3' x4' x5']; % data must be organized column-wise

u=y(nDiscard+1:N,:);

end

