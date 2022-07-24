%% ARFITCAPS 
%        Capsule function to call the arfit.m routine, part of 
%        "ARfit: Multivariate Autoregressive Model Fitting" package.
%
%% Syntax:
%      [pf,A,ef] = ARFITCAPS(u,IP)
%
%% Input Arguments
%   u:          time series
%   IP:         VAR model order
%
%% Output Arguments 
%   pf:         covariance matrix provided by ARFIT routine
%   A:          AR estimate matrix provided by ARFIT routine
%   ef:         forward residuals provided by ARRES routine
%
%% Description:
% ARFITCAPS is capsule that calls arfit.m and arres.m routines, part of 
% Autoregressive Model Fitting" package, which implements algorithms 
% as described in the following articles:
%
%    [1] Neumaier A & Schneider T, 2001. Estimation of parameters and
%         eigenmodes of multivariate autoregressive models. ACM Trans Math
%         Softw, 27:27-57.
%    [2] Schneider T & Neumaier A, 2001. Algorithm 808: ARfit - A Matlab
%        package for the estimation of parameters and eigenmodes of multivariate
%        autoregressive models. ACM Trans Math Softw, 27:58-65. 
%
%  If you are interested in using ARfit algorithm for VAR model estimation, we
%  advise you to get the software from Tapio Schneider's website at
%                 https://climate-dynamics.org/software/#arfit
%  or from Mathworks.com File Exchange site (verified on August 27, 2021, but
%  not tested)
%        https://www.mathworks.com/matlabcentral/fileexchange/174-arfit,
%
%  and, before using it, verify the license terms, it seems to be a copyrighted
%  material by the Association for Computing Machinery, Inc.
%
%%  Note: As described by the authors, acf.m in ARfit needs Signal Processing 
%  Toolbox (TM), as it requires XCORR, a cross-correlation function estimator.
%
%  ARfit availability was checked on August 13, 2015, and August 27, 2021. KS
%
%  The version we have tested and included was obtained on February 24, 2011 from
%       www.gps.caltech.edu/~tapio/arfit/index.html, 
%  which is now obsolete.
%
%% See also ARFIT, MVAR, MCARNS, MCARVM, CMLSM
%          arfit | <arfitcaps.html> |<mvar.html>|

% (C) Koichi Sameshima & Luiz A. Baccala, 2021. 
% See file license.txt in installation directory for licensing terms.

%%

function [pf,A,ef] = arfitcaps(u,IP)

if ~exist('arfit.m','file')
   help arfitcaps
   error('ARfit.m not found. Get the ARfit package from Tapio Schneider''s web site.')
end;

v = u';
[w, Au, C, sbc, fpe, th] = arfit(v,IP,IP);
pf = C;

if IP >= 20
   [siglev,res] = arres(w,Au,v,IP+1);
else
   [siglev,res] = arres(w,Au,v);
end;

% Variable 'siglev' is not used.

ef = res'; 
A = zeros(length(w),length(w),IP);
for i = 1:IP
   A(:,:,i) = Au(:,(i-1)*length(w)+1:i*length(w));
   wu = ceil(length(ef)*rand(size(w)));
   if length(ef)<length(v)
      ef = [ef ef(:,wu(1))];
   else
      ef = ef(:,1:length(v));
   end
end
