%% Daniel Y. Takahashi  - Medidas de Fluxo de Informacao com Aplicacao em Neurociencia (Thesis in Portuguese)
%   [http://www.princeton.edu/~dtakahas/publications/TeseDYTFinal.pdf]
%
% Thesis example: An invertible VMA process with an infinite order VAR process
% representation was simulated

% (C) Daniel Y. Takahashi, Koichi Sameshima & Luiz A. BaccalÁ, 2022. 
% See file license.txt in installation directory for licensing terms.

function u = takahashi_thesis_dat(nPoints)
N=30000;
disp('======================================================================');
disp('               VAR[2] Trivariate Model (Example 8.1.1, page 153)');
disp('                        Thesis D.Y.Takahashi (2008)')
disp('                             Y --> X');
disp('======================================================================');
randn('state', sum(100*clock));
w1 = randn(1,N);
w2 = randn(1,N);
w3 = randn(1,N);

X = zeros(1,N);
Y = zeros(1,N);
Z = zeros(1,N);

for k = 1:2
   X(k) = w1(k);
   Y(k) = w2(k);
   Z(k) = w3(k);
end
% VAR model 
for k = 3:30000
   X(k) = -0.64*Y(k-2) + 0.8*Z(k-1) + w1(k);
   Y(k) = w2(k);
   Z(k) = 0.8*Y(k-1) + w3(k);
end

yy = [X' Y' Z']; % data must be organized column-wise
nDiscard = 1000; % number of points discarded at beginning of simulated series
if nargin < 1
  nPoints = 200;  % number of analyzed samples points
end

u = yy(nDiscard+1:nDiscard+nPoints,:);

disp('======================================================================');
disp('                      End of Takahashi VAR model');
disp('======================================================================');

