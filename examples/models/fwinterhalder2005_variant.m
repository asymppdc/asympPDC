function [ u ] = fwinterhalder2005_variant( nPoints, nDiscard )

% Winterhalder et al. Comparison of linear signal processing techniques to
% infer directed interactions in multivariate neural systems. 
% Signal Processing 85:2137--2160, 2005.
%    [http://dx.doi.org/10.1016/j.sigpro.2005.07.011]
%
% Example: Seven random independent variables 


disp('======================================================================');
disp('       Winterhalder et al. Signal Processing. 85:2137--60, 2005')
disp('        Variant of Random Independent Process with 7 variables')
disp(' [sigma1=500|sigma2=1|sigma3=500|sigma4=1|sigma5=1|sigma6=1|sigma7=1]');
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

N = nDiscard+nPoints; % number of simulated points

randn('state', sum(100*clock))
ei=randn(10,N);

sigma1=500; sigma2=1; sigma3=500; sigma4=1; sigma5=1; sigma6=1; sigma7=1;

x1=sigma1*ei(1,:); x2=sigma2*ei(4,:); x3=sigma3*ei(6,:); x4=sigma4*ei(7,:);
x5=sigma5*ei(5,:); x6=sigma6*ei(2,:); x7=sigma7*ei(8,:);

y=[x1' x2' x3' x4' x5' x6' x7'];

u=y(nDiscard+1:N,:);

end

