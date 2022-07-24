function [ u ] = fguo2008_linear( nPoints, nDiscard, flgManual )

% Guo, Wu, Ding & Feng. Uncovering interactions in frequency domains.
%         PLoS Computational Biology, 4(5):1-10, February 8, 2008.
%           [http://dx.plos.org/10.1371/journal.pcbi.1000087] 
% 
% Page 2 Toy Model Example 5-dimensional VAR[3] with large common exogenous input

disp(repmat('=',1,100))
disp('                   Five dimensional linear VAR[3] Model')
disp('                    with large common exogenous input')
disp('                      Guo et al. February 8, 2008.')
disp('                 x1==>x2  x1-->x3 x1==>x4 x4-->x5 x5-->x4 ');
disp('             Instantaneous causality in between all variables.');
disp(repmat('=',1,100))

%flgManual = 0;

if (nargin == 0),
   nPoints = 1000;
   nDiscard = 5000;
   disp(['Adopting default ' int2str(nDiscard) ' discarding points, and'])
   disp(['generating ' int2str(nPoints) ' simulation data point.'])
   
elseif   (nargin < 2),
   nDiscard = 5000;
   disp(['Adopting default ' int2str(nDiscard) ' discarding points.'])
elseif (nargin < 3)
   flgManual = 0;
end;

if (nDiscard < 1),
   nDiscard = 5000;
   disp(['Adopting default ' int2str(nDiscard) ' discarding points.'])
end;

if nPoints < 10,
   nPoints = 100;
   disp(['Adopting default ' int2str(nPoints) ' simulation data points.'])
end;
%
randn('state', sum(100*clock))

%randn('state', 100*pi); % This command line was used to repeat 
%                        % the analysis with the same data set to check 
%                        % alg and criterion selection.

%% Variables initialization

N = nDiscard+nPoints; % number of simulated points

if ~exist('flgManual','var')
   flgManual = 0
end;

if flgManual
   if exist('randn_manual_state.mat','file'),
      load randn_manual_state
      randn('state', s);
      disp(['Using state saved in "randn_manual_state.mat" to reproduce figure in the manual.']);
   else
      disp(['Did not find "randn_manual_state.mat" file.' ...
         'Assigned "sum(100*clock)" initial state.']);
      randn('state', sum(100*clock))
   end;
else
   randn('state', sum(100*clock))
   disp('Assigned "sum(100*clock)" initial state.2')
end;

% Variables initialization
ei = randn(7,N);
x1 = zeros(1,N);
x2 = zeros(1,N);
x3 = zeros(1,N);
x4 = zeros(1,N);
x5 = zeros(1,N);

a1 = rand; a2 = rand; a3 = rand; a4 = rand; a5 = rand;
b1 = 2; b2 = 2; b3 = 2; b4 = 2; b5 = 2;
c1 = 5; c2 = 5; c3 = 5; c4 = 5; c5 = 5;
%%
% The model with all above parameters set to zero corresponds to Example 4 
% of Baccala & Sameshima (2001). See example baccala2001a_ex4.m
% a1 = 0; a2 = 0; a3 = 0; a4 = 0; a5 = 0;
% b1 = 0; b2 = 0; b3 = 0; b4 = 0; b5 = 0;
% c1 = 0; c2 = 0; c3 = 0; c4 = 0; c5 = 0;

%%
%

disp('Model parameters:')
disp(sprintf('a1=%f; a2=%f; a3=%f; a4=%f; a5=%f;',a1,a2,a3,a4,a5))
disp('b1 = b2 = b3 = b4 = b5 = 2;')
disp('c1 = c2 = c3 = c4 = c5 = 5;')
disp(' ')
disp(repmat('=',1,100))

for t = 1:4,
   x1(t) = randn(1); x2(t) = randn(1); x3(t) = randn(1); x4(t) = randn(1);
   x5(t) = randn(1); 
end;

chLabels = {'x_1';'x_2';'x_3';'x_4';'x_5'}; %or %chLabels = []; 

for t = 5:N,
   x1(t) = 0.95*sqrt(2)*x1(t-1) - 0.9025*x1(t-2) + ei(1,t) + ...
                                 a1*ei(6,t) + b1*ei(7,t-1) + c1*ei(7,t-2);
   x2(t) = 0.5*x1(t-2) + ei(2,t) + ...
                                 a2*ei(6,t) + b2*ei(7,t-1) + c2*ei(7,t-2);
   x3(t) =-0.4*x1(t-3) + ei(3,t) + ...
                                 a3*ei(6,t) + b3*ei(7,t-1) + c3*ei(7,t-2);
   x4(t) =-0.5*x1(t-2) + 0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) ...
                     + ei(4,t) + a4*ei(6,t) + b4*ei(7,t-1) + c4*ei(7,t-2);
   x5(t) =-0.25*sqrt(2)*x4(t-1) + 0.25*sqrt(2)*x5(t-1) ...
                     + ei(5,t) + a5*ei(6,t) + b5*ei(7,t-1) + c5*ei(7,t-2);
end;

y = [x1' x2' x3' x4' x5']; % data must be organized row-wise
u = y(nDiscard+1:N,:);

[nSegLength,nChannels] = size(u);

end
