%% ASYMP_PDC
%        Compute complex and squared partial directed coherence for any of three
%        metrics -- Euclidean, diagonal and informational -- along with
%        asymptotic statistics from vector autoregressive (VAR) coefficients in
%        the frequency domain.
%
%% Syntax
%        c = ASYMP_PDC(u,A,pf,nFreqs,metric,alpha)
%
%% Input Arguments
%        u:          multiple row vectors time series
%        A:          AR estimate matrix by MVAR routine
%        pf:         covariance matrix provided by MVAR routine
%        nFreqs:     number of point in [0,fs/2) frequency scale range
%        metric:     'euc'  -- Euclidean   ==> original PDC
%                    'diag' -- diagonal    ==> gPDC (generalized)
%                    'info' -- information ==> iPDC
%        alpha:      significance level
%                    if alpha is zero, statistical analysis will not be
%                    performed for faster computation.
%
%% Output Arguments 
%        c: struct variable with following fields:
%        |-- .pdc       - PDC complex estimates
%        |-- .pdc2      - |PDC|^2 estimates
%        |-- .pvalues   - p-values associated to PDC2 estimates. 
%        |-- .th        - |PDC|^2 threshold value with (1-avalue) significance 
%        |                level.
%        |-- .{ci1,ci2} - lower and upper (1 - alpha) confidence interval of 
%        |                |PDC|^2 estimates
%        |-- .metric    - metric used for PDC calculation 
%        |-- .alpha     - significance level
%        |-- .p         - VAR model order
%        |-- .patdenr   - 
%        |-- .patdfr    - degree of freedom 
%        |-- .SS        - power spectra
%        +-- .coh2      - squared spectral coherence, 
%      i.e.
%        .{pdc,pdc2,pvalues,th,ci1,ci2,metric,alpha,p,patdenr,patdfr,SS,coh2}
%
%% Description
%   Compute all three types of PDC --- connectivity measure --- and allied
%   asymptotic statistics [2] measures for the chosen metric option: 
%        * 'euc'  - original or Euclidean PDC as proposed in [1]; 
%        * 'diag' - generalized PDC; 
%        * 'info' - information PDC. 
%
%% Example:
% 
% Annual sunspot numbers  and the melanoma cases (x10^5) in the State of
% Connecticuts, USA, from 1936 to 1972, given by
%
%   % -------- Example script start here---------------------------------------
%
%   % Data and pre-processing
%   u = [ 40 115 100  80  60  40  23  10  10  25  75 145 130 130  80  65  20 ...
%         10   5  10  60 190 180 175 120  50  35  20  10  15  30  60 105 105 ...
%         105  80  65; ...
%        0.9 0.8 0.8 1.3 1.4 1.2 1.7 1.8 1.6 1.5 1.5 2.0 2.5 2.7 2.9 2.5 3.1 ...
%        2.4 2.2 2.9 2.5 2.6 3.2 3.8 4.2 3.9 3.7 3.3 3.7 3.9 4.1 3.8 4.7 4.4 ...
%        4.8 4.8 4.8];
%   [nChannels,nSegLength] =size(u);
%   if nChannels > nSegLength, u = u.'; [nChannels,nSegLength]=size(u); end;
%   for i=1:nChannels, u(i,:)=detrend(u(i,:)); end;
%
%   % Pre-calculated AR model coefficients, A and pf
%   A        = [ 0.8280  -12.1097; 0.0009   -0.1258]; % VAR model estimate.
%   A(:,:,2) = [-0.0724   13.4798; 0.0036   -0.1391];
%   A(:,:,3) = [-0.3561  -36.4805; 0.0013   -0.0735];
%   pf = [568.0873  -1.5815; -1.5815 0.0474];
%
%   % Analysis parameters
%   nFreqs = 128;    % number of points to calculate PDC in the frequency scale
%   measure = 'pdc'; % Connectivity measure
%   metric = 'info'; % calculating information PDC
%   alpha  = 0.01;   % significance level
%   c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Calculate |iPDC|^2 with 
%                                              % statistics
%
%   % Let's visualize PDC plots with default parameters. Try this first.
%   xplot([],c);  
%
%   % -------- Example end here------------------------------------------------
%   %
%   % -------- Further xplot options example start here------------------------
%
%   % Then let's xplot PDC2 with threshold, confidence interval, and labels, 
%   % also with spectral coherence. 
%   % First, let's set some "xplotting" parameters:
%
%   chLabels    = {'Sunspot';'Melanoma'}; % Channel labels
%   flgColor    = 0;      fs = 1;    w_max = 0.5; 
%   flgPrinting = [1 1 1 2 2 1 1]; flgScale = 1; 
%   flgMax      = 'all';    flgSignifColor = 4;
%   xplot('THIS IS A WINDOW TITLE',c,flgPrinting,fs,fs/2,chLabels,flgColor, ...
%                                               flgScale,flgMax,flgSignifColor);
%   % and let's add an informative title
%   xplot_title(alpha,metric,measure,'Sunspot-Melanoma Example')
%
%   % -------- Second part of xplot example end here---------------------------
%
%% References
%
% [1] L.A. Baccala and K. Sameshima. Partial directed coherence: a new concept
%     in neural structure determination. Biol Cybern 84:463--474,2001.
%     <https://doi.org/10.1007/PL00007990>
%
% [2] D.Y. Takahashi, L.A.B. Baccala and K. Sameshima, Connectivity inference
%     between neural structures via partial directed coherence. J Appl Stat
%     34:1259--1273, 2007. <https://doi.org/10.1080/02664760701593065>
%
% [3] L.A. Baccala, C.S.N. De Brito, D.Y. Takahashi and K. Sameshima. Unified
%     asymptotic theory for all partial directed coherence forms. Philos T Roy
%     Soc A 371(1997):1--13, 2013. <https://doi.org/10.1098/rsta.2012.0158>
%
% See also PDC_ALG, ASYMP_DTF, MVAR, MCARNS, MCARVM, CMLSM, ARFIT
%          | <asymp_pdc.html> |<asymp_dtf.html>|

% (C) Koichi Sameshima & Luiz A. Baccala, 2021. See file license.txt in
%     the installation directory for licensing terms.

function c = asymp_pdc(u,A,pf,nFreqs,metric,alpha)

if nargin < 6
    error('ASYMP_PDC requires six input arguments.')
end
[m,n] = size(u);
if m > n
    u = u.';
end
np = length(u);
[nChannels, ~, p] = size(A);  % ~ = dummy variable == nChannels
Af = A_to_f(A, nFreqs);

flgVerbose = 0;

% Variables pre-alocation
pdc = zeros(nChannels,nChannels,nFreqs);  % To hold complex PDC
pdc2 = zeros(nChannels,nChannels,nFreqs); % for squared-PDC

if alpha ~= 0
   th      = zeros(nChannels,nChannels,nFreqs);
   ci1     = zeros(nChannels,nChannels,nFreqs);
   ci2     = zeros(nChannels,nChannels,nFreqs);
   varass1 = zeros(nChannels,nChannels,nFreqs);
   varass2 = zeros(nChannels,nChannels,nFreqs);
   patdfr  = zeros(nChannels,nChannels,nFreqs);
   patdenr = zeros(nChannels,nChannels,nFreqs);
   pvalues = zeros(nChannels,nChannels,nFreqs);
   if flgVerbose
      switch lower(metric)
         case {'euc'}
            disp('* Original PDC and asymptotic statistics')
         case {'diag'}
            disp('* Generalized PDC and asymptotic statistics')
         case {'info'}
            disp('* Information PDC and asymptotic statistics')
         otherwise
            error('Unknown metric.')
      end
   end
elseif flgVerbose
   switch lower(metric)
      case {'euc'}
         disp('* Original PDC estimation')
      case {'diag'}
         disp('* Generalized PDC estimation')
      case {'info'}
         disp('* Information PDC estimation')
      otherwise
         error('Unknown metric.')
   end
end

% ======================================================================

switch lower(metric)
   case {'euc'}           % for original PDC
      dpdc_dev = zeros(1,(nChannels*(nChannels + 1))/2);
      pfe = eye(nChannels);  % for complex PDC calculation  
      
   case {'diag'}          % for gPDC
      evar_d = mdiag(pf);
      evar_d_big = kron(eye(2*nChannels), evar_d);
      pinv_evar_d_big = pinv(evar_d_big);

      %'derivative of vec(Ed-1) by vecE'
      de_deh = Dup(nChannels);
      debig_de = fdebig_de(nChannels);
      
      dedinv_dev = diagtom(vec(-pinv_evar_d_big*pinv_evar_d_big));
      dedinv_deh = dedinv_dev*debig_de*de_deh;
      pfe = pf;          % for complex gPDC calculation
      
   case {'info'}         % for iPDC
      evar_d = mdiag(pf);
      evar_d_big = kron(eye(2*nChannels), evar_d);
      pinv_evar_d_big = pinv(evar_d_big);
      
      evar_big = kron(eye(2*nChannels), pf);
      pinv_evar_big = sparse(pinv(evar_big));

      %'derivative of vec(Ed-1) by vecE'
      de_deh = Dup(nChannels);
      debig_de = fdebig_de(nChannels);
      
      dedinv_devd = sparse(diagtom(vec(-pinv_evar_d_big*pinv_evar_d_big)));
      dedinv_dehd = sparse(dedinv_devd*debig_de*de_deh);
      
      %dedinv_dev = sparse(-kron(inv_e.', inv_e));
      dedinv_dev = sparse(-kron(pinv_evar_big.', pinv_evar_big));
      dedinv_deh = sparse(dedinv_dev*debig_de*de_deh);
      pfe = pf;   % necessary for complex iPDC calculation

   otherwise
      error('Unknown metric.')
end

gamma = bigautocorr(u, p);
omega = kron(inv(gamma), pf);
omega_evar = 2*pinv(Dup(nChannels))*kron(pf, pf) * pinv(Dup(nChannels)).';

%if (isOctave())
icdf_norm_alpha = norminv(1 - alpha/2.0,0,1);
%else
%   icdf_norm_alpha = icdf('norm',1 - alpha/2.0,0,1);
%end

for ff = 1:nFreqs
   f = (ff - 1)/(2*nFreqs);     %Corrected 7/25/2011, f starting at 0 rad/s.
   Ca = fCa(f, p, nChannels);
   
   a = Af(ff,:,:); a = a(:);
   a = [real(a); imag(a)];

   omega2 = Ca*omega*Ca.';
   L = fChol(omega2);

   for j = 1:nChannels
      Ij = fIj(j,nChannels);
      switch lower(metric)
         case {'euc'}           % for PDC
            Ije = Ij;
            
         case {'diag'}          % for gPDC
            Ije = Ij*pinv_evar_d_big;
            
         case {'info'}          % for iPDC
            Ije = Ij*pinv_evar_big*Ij;
            
         otherwise
            error('Unknown metric.')
      end
      
      for i = 1:nChannels
         Iij = fIij(i,j,nChannels);
         %For metric = diag or info case, include evar in the expression'
         switch lower(metric)
            case {'euc'}           % PDC
               Iije = Iij;
               
            case {'diag'}          % gPDC
               Iije = Iij*pinv_evar_d_big;
               
            case {'info'}          % iPDC
               Iije = Iij*pinv_evar_d_big;
               
            otherwise
               error('Unknown metric.')
         end

         num = a.'*Iije*a;
         den = a.'*Ije*a;
         pdc2(i,j,ff) = num/den; % $|PDC_{ij}(ff)|^2$
         
         numd = Af(ff,i,j)/sqrt(pfe(i,i));
         pdc(i,j,ff) = numd/sqrt(den); % complex PDC calculation
         
         % If alpha == 0, no statistical calculation for faster computation.
         if alpha ~= 0
            %'Add evar derivation'
            switch lower(metric)
               case {'euc'}
                  %nop
                  
               case {'diag'}
                  %'derivative of num by vecE'
                  dnum_dev = kron((Iij*a).', a.')*dedinv_deh;
                  %'derivative by den by vecE'
                  dden_dev = kron((Ij*a).', a.')*dedinv_deh;
                  dpdc_dev = (den*dnum_dev - num*dden_dev)/(den^2);
                  
               case {'info'}
                  %'derivative of num by vecE'
                  dnum_dev = kron((Iij*a).', a.')*dedinv_dehd;
                  %'derivative of den by vecE'
                  dden_dev = kron((Ij*a).', a.'*Ij)*dedinv_deh;
                  dpdc_dev = (den*dnum_dev - num*dden_dev)/(den^2);
                  
               otherwise
                  error('Unknown metric.')
            end

            G1a = 2 * a.' * Iije/den - 2*num*a.'*Ije/(den^2);
            G1 = -G1a * Ca;
            varalpha = G1 * omega * G1.';
            varevar = dpdc_dev * omega_evar * dpdc_dev.';
            varass1(i,j,ff) = (varalpha + varevar)/np;
            
            ci1(i,j,ff) = pdc2(i,j,ff) ...
                    - sqrt(varass1(i,j,ff))*icdf_norm_alpha;
            ci2(i,j,ff) = pdc2(i,j,ff) ...
                    + sqrt(varass1(i,j,ff))*icdf_norm_alpha;

            G2a = Iije/den;
                        
            d = fEig(real(L), real(G2a)); % real() 28May2013            

            patdf = (sum(d).^2)./sum(d.^2);
            patden = sum(d)./sum(d.^2);

            th(i,j,ff) = chi2inv((1 - alpha), patdf)/(patden*np);
            
            pvalues(i,j,ff) = 1 - chi2cdf(pdc2(i,j,ff)*patden*np, patdf);
 
            varass2(i,j,ff) = patdf/(patden*np).^2;
            patdfr(i,j,ff) = patdf;
            patdenr(i,j,ff) = patden;
         else % alpha is zero, no asymptotics calculation
            % nop;
         end
      end  % j = nChannels
   end  % i = nChannels
end  % ff = freqs

% assigning values to the c structure.
if alpha ~= 0
   c.pdc = pdc;   % Complex PDC2/gPDC2/iPDC2
   c.pdc2 = pdc2;
   c.th = th;
   c.ci1 = ci1;
   c.ci2 = ci2;
   c.metric = metric;
   c.alpha = alpha;
   c.p = p;
   c.pvalues = pvalues; % p-values associated to PDC2/gPDC2/iPDC2
   c.patden = patdenr;
   c.patdf = patdfr;
   c.varass1 = varass1;
   c.varass2 = varass2;   
%  Statistically significant PDC2 on frequency scale
   pdc2_temp = ((abs(pdc2) - abs(th)) > 0).*pdc2 ...
                                           + ((abs(pdc2) - abs(th)) <= 0)*(-1);
   pdc2_temp(pdc2_temp < 0) = NaN;
   c.pdc2_th = pdc2_temp;
   
else % No statistics, just the |PDC^2| values/
   c.pdc = pdc;
   c.pdc2 = pdc2;
   c.metric = metric;
   c.alpha = 0;
   c.p = p;
   c.pvalues = [];   
   c.th = [];
   c.ci1 = [];
   c.ci2 = [];
   c.patdenr = [];
   c.patdfr = [];
   c.varass1 = [];
   c.varass2 = [];
end

% Power spectra and coherence assignments to c struct.
c.SS = ss_alg(A, pf, nFreqs);
c.coh2 = coh_alg(c.SS);
%
end

%% 
% Subfunctions
%

%==========================================================================
function gamma = bigautocorr(x, p)
% Autocorrelation. Data in row-wise orientation. From order 0 to p-1.
% Output: n x n blocks of autocorr of lags i. (Nuttall Strand matrix)'''

[n, nd] = size(x);

gamma = zeros(n*p, n*p);
for i = 1:p
   for j = 1:p
      gamma(((i - 1)*n + 1):i*n, ((j - 1)*n + 1):j*n) = ...
                                           xlag(x, i - 1)*(xlag(x,j - 1).')/nd;
   end
end
end

%==========================================================================
function c = xlag(x,tlag)
if tlag == 0
   c = x;
else
   c = zeros(size(x));
   c(:,(tlag + 1):end) = x(:,1:(end - tlag));
end
end

%==========================================================================
function d = fEig(L, G2)
%'''Returns the eigenvalues'''

%L = mat(cholesky(omega, lower=1))
D = L.'*G2*L;

%disp('fEig: eig or svd?')

d  = svd(D);
d1 = sort(d);
%
% the two biggest eigenvalues no matter which values (non negative by
% construction
%
d = d1(length(d) - 1:length(d));

if (size(d) > 2)
   disp('More than two Chi-squares in the sum:')
end
end

%==========================================================================
function c = fIij(i,j,n)
%'''Returns Iij of the formula'''
Iij = zeros(1,n^2);
Iij(n*(j - 1) + i) = 1;
Iij = diag(Iij);
c = kron(eye(2), Iij);

c = sparse(c);  % SPARSED
end

%==========================================================================
function c = fIj(j,n)
%'''Returns Ij of the formula'''
Ij = zeros(1,n);
Ij(j) = 1;
Ij = diag(Ij);
Ij = kron(Ij,eye(n));
c = kron(eye(2), Ij);

c = sparse(c);  % SPARSED
end

%==========================================================================
function d = fCa(f, p, n)
%'''Returns C* of the formula'''
C1 = cos(-2*pi*f*(1:p));
S1 = sin(-2*pi*f*(1:p));
C2 = [C1; S1];
d = kron(C2, eye(n^2));
end

%==========================================================================
function c = fdebig_de(n)
%'''Derivative of kron(I(2n), A) by A'''
%c = kron(TT(2*n, n), eye(n*2*n)) * kron(eye(n), kron(vec(eye(2*n)), eye(n)));
A = sparse(kron(TT(2*n, n), eye(n*2*n)));
B = sparse(kron(vec(eye(2*n)), eye(n)));
c = A * kron(eye(n), B);
c = sparse(c);
end

%==========================================================================
function c = vec(x)
%vec = lambda x: mat(x.ravel('F')).T
c = x(:);
end

%==========================================================================
function t = TT(a,b)
%''' TT(a,b)*vec(B) = vec(B.T), where B is (a x b).'''
t = zeros(a*b);
for i = 1:a
   for j = 1:b
      t((i - 1)*b + j,(j - 1)*a + i) = 1;
   end
end
t = sparse(t);
end

%==========================================================================
function L = fChol(omega)
% Try Cholesky factorization
try
   L = chol(omega)';
   % If there's a small negative eigenvalue, diagonalize
catch
   %   disp('linalgerror, probably IP = 1.')
   [v,d] = eig(omega);
   L = zeros(size(v));
   for i = 1:length(d)
      if d(i,i) < 0
         d(i,i) = eps;
      end
      L(:,i) = v(:,i)*sqrt(d(i,i));
   end
end
end

%==========================================================================
function c = diagtom(a)
b = sparse(a');
c = sparse(diag(b(:)));
end

%==========================================================================
function c = mdiag(a)
%  diagonal matrix
c = diag(diag(a)); % UNSPARSED
end

%==========================================================================
function d=Dup(n)
% '''D*vech(A) = vec(A), with symmetric A'''
d = zeros(n*n, (n*(n + 1))/2);
count = 1;
for j = 1:n
   for i = 1:n
      if i >= j
         d((j - 1)*n + i,count) = 1;
         count = count + 1;
      else
         d((j - 1)*n + i,:) = d((i - 1)*n + j,:);
      end
   end
end
end

%==========================================================================

%%
%
%        1         2         3         4         5         6         7         8         9
%23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890

%% Change Log:
%
% [20110725]: The asymp_pdc routine was corrected to match the frequency 
%             range with plotting routine, f = 0 was included in the "frequency" 
%              for-loop:
%                                for ff = 1:nFreqs,
%                                   f = (ff-1)/(2*nFreqs); %
%                                        ^?^^
% [20150107]: Optimization \(^o^)/

% [EOF]