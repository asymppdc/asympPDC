function [ u aState] = fschelter2006( nPoints, nDiscard, flgRepeat)


% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.

if nargin < 3, flgRepeat = 0; end; 

disp('======================================================================');
disp('         Schelter et al. J Physiology - Paris 99:37-46, 2006.')
disp('                Linear penta-variate VAR[4]-process')
disp('       x1-->x2  x1-->x4 x2-->x4 x4==>x5 x5-->x1  x5-->x2 x5-->x3 ');
disp('======================================================================');

if flgRepeat,
   if exist('schelter2006_state.mat', 'file') == 2,
      load schelter2006_state.mat % Retrieving saved state number.
      if ~exist('aState','var'),
         disp('State number does not exist. State number initialized.')
         aState = sum(100*clock);
      else
         disp('Using saved state number to repeat simulation with same dataset.');
      end;
   else
      disp('File schelter2006_state.mat not found. Using newly state number.')
      aState = sum(100*clock);
      save schelter2006_state.mat aState
   end
else
   aState = sum(100*clock); % Starting a new state number.
   
end;


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


N = nDiscard + nPoints; % number of simulated points

randn('state',aState)

ei=randn(5,N);
x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);
x4=zeros(1,N);
x5=zeros(1,N);

% Variables initialization
for t=1:6,
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1); x4(t)=randn(1);
   x5(t)=randn(1);
end;

for t=7:N,
   x1(t) = 0.4*x1(t-1) - 0.5*x1(t-2) + 0.4*x5(t-1) + ei(1,t);
   x2(t) = 0.4*x2(t-1) - 0.3*x1(t-4) + 0.4*x5(t-2) + ei(2,t);
   x3(t) = 0.5*x3(t-1) - 0.7*x3(t-2) - 0.3*x5(t-3) + ei(3,t);
   x4(t) = 0.8*x4(t-3) + 0.4*x1(t-2) + 0.3*x2(t-2) + ei(4,t);
   x5(t) = 0.7*x5(t-1) - 0.5*x5(t-2) - 0.4*x4(t-1) + ei(5,t);
end;

y = [x1' x2' x3' x4' x5']; % data must be organized row-wise
u = y(nDiscard+1:N,:);

end

