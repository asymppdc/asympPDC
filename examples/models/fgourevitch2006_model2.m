function [ u ] = fschelter2006_model2( nPoints, nBurnIn)

% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.

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


N = nBurnIn + nPoints; % number of simulated points

disp(repmat('=',1,100))
disp('             Gourevitch et al. Biol Cybern 95:349-69, 2006.')
disp('     Model 2: Linear Bivariate model with Bidirectional Influence ')
disp('    at different frequencies in each direction, with Common Source')
disp('         x1-->x2    x2-->x1   x1--x2 (With instantaneous causality)');
disp(repmat('=',1,100))

randn('state', sum(100*clock))
wi=randn(3,N); % ws(t) = wi(3,t)
x1=zeros(1,N);
x2=zeros(1,N);
for t=1:3
   x1(t)=wi(1,t);
   x2(t)=wi(2,t);
end
for t=4:N
   x1(t) =  0.95*sqrt(2)*x1(t-1) - 0.9025*x1(t-2) ... 
                                 - 0.9*x2(t-1) + 0.5*wi(1,t) + 0.5*wi(3,t);
   x2(t) = -1.05*x2(t-1) - 0.85*x2(t-2) ...
                                 - 0.8*x1(t-1) + 0.5*wi(2,t) + 0.5*wi(3,t);
end

y=[x1' x2'];     % data must be organized column-wise
u = y(nBurnIn+1:N,:);

end

