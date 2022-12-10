%% ASYMP_DTF
%        Compute DTF connectivity measures magnitude, from series j--> i, for
%        any of three of metrics -- Euclidean, diagonal and information --
%        as well as asymptotic statistics from vector autoregressive (VAR)
%        coefficients in the frequency domain.
%
%% Syntax:
%        c = ASYMP_DTF(u,A,pf,nFreqs,metric,alpha)
%
%% Input Arguments:
%        u      - multiple row vectors time series
%        A      - AR estimate matrix obtained via MVAR routine
%        pf     - covariance matrix provided via MVAR routine
%        nFreqs - number of point in [0,fs/2) frequency scale
%        metric - 'euc':  -- Euclidean   ==> original DTF
%                 'diag': -- diagonal    ==> DC (directed coherence)
%                 'info': -- information ==> iDTF
%        alpha  - significance level
%                 if alpha is zero, statistical analysis won't be performed
%
%% Output Arguments:
%        c struct variable with following fields:
%        |-- .dtf       - complex DTF estimates
%        |-- .dtf2      - |DTF|^2 estimates
%        |-- .pvalues   - p-values associated to DTF2 estimates.
%        |-- .th        - |DTF|^2 threshold value with (1-alpha) significance
%        |                level.
%        |-- .{ci1,ci2} - lower and upper (1 - alpha) confidence interval of
%        |                |DTF|^2 estimates
%        |-- .metric    - metric used for DTF calculation
%        |-- .alpha     - significance level
%        |-- .p         - VAR model order
%        |-- .patdenr   -
%        |-- .patdfr    - degree of freedom
%        |-- .SS        - power spectra
%        +-- .coh2      - squared spectral coherence
%    or
%        c.{dtf,dtf2,pvalues,th,ci1,ci2,metric,alpha,p,patdenr,patdfr,SS,coh2}
%
%% Description:
%   Compute all three types of DTF --- Granger influentiability measure and
%   their allied statistical measures of asymptotic statistics for metric
%   option:
%        * 'euc'  - original or Euclidean DTF as proposed in [Kaminski &
%                   Blinowska, 2001];
%        * 'diag' - Directed Coherence (DC) or gDTF (generalized);
%        * 'info' - information DTF.
%
%% Example:
%
% Annual sunspot numbers  and the melanoma cases (10^5) in the State of
% Connecticuts, USA, from 1936 to 1972, given by
%
%   u = [ 40 115 100  80  60  40  23  10  10  25  75 145 130 130  80  65  20 ...
%         10   5  10  60 190 180 175 120  50  35  20  10  15  30  60 105 105 ...
%         105  80  65; ...
%        0.9 0.8 0.8 1.3 1.4 1.2 1.7 1.8 1.6 1.5 1.5 2.0 2.5 2.7 2.9 2.5 3.1 ...
%        2.4 2.2 2.9 2.5 2.6 3.2 3.8 4.2 3.9 3.7 3.3 3.7 3.9 4.1 3.8 4.7 4.4 ...
%        4.8 4.8 4.8];
%
%
%   [nChannels,nSegLength] =size(u);
%   if nChannels > nSegLength, u = u.'; [nChannels,nSegLength]=size(u); end;
%   for i=1:nChannels, u(i,:)=detrend(u(i,:)); end;
%
%   A        = [ 0.8280  -12.1097; 0.0009   -0.1258]; % VAR model estimate.
%   A(:,:,2) = [-0.0724   13.4798; 0.0036   -0.1391];
%   A(:,:,3) = [-0.3561  -36.4805; 0.0013   -0.0735];
%
%   pf = [568.0873  -1.5815; -1.5815 0.0474];
%
%   nFreqs = 128;     % number of points to calculate DTF in the frequency scale
%   metric = 'info';  % calculating information DTF
%   alpha  = 0.01;    % significance level
%   c = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Calculate iDTF with statistics
%
%   chLabels={'Sunspot';'Melanoma'}; % Channel labels
%   flgColor    = 0;      fs = 1;    w_max = 0.5;
%   flgPrinting =[1 1 1 2 2 0 1]; flgScale = 1;
%   flgMax      = 'all';    flgSignifColor = 3;
%   figure; xplot(c) % Visualiza PDC plots. Try this first.
%
%   figure;
%   xplot(c,flgPrinting,fs,w_max,chLabels,flgColor); % DTF plot with confidence
%                                                    % interval
%   figure;
%   xplot(c,flgPrinting,fs,fs/2,chLabels,flgColor,flgScale,flgMax, ...
%                                                              flgSignifColor);
%
%% References:
%   [1] M.J. Kaminski and K.J. Blinowska. A new method of the description of the
%   information flow in the brain structures. Biol Cybern 65:203--210,1991.
%   <https://doi.org/10.1007/bf00198091>
%
%   [2] L.A.B. Baccala, D.Y. Takahashi and K. Sameshima. Directed transfer
%   function: unified asymptotic theory and some of its implications. IEEE T
%   Bio-Med Eng 63:2450--2460, 2016.
%   <https://doi.org/10.1109/TBME.2016.2550199>
%
%   [3] F. Rezaei, O. Alamoudi, S. Davani and S. Hou  (2022). Fast
%   Asymptotic Algorithm for Real-Time Causal Connectivity Analysis of
%   Multivariate Systems and Signals. Signal Processing.
%   <https://doi.org/10.1016/j.sigpro.2022.108822>
%
%% See also: DTF_ALG, ASYMP_PDC, MVAR, MCARNS, MCARVM, CMLSM, ARFIT, FASTASYMPALG, FASTASYMPALG0, FASTASYMPALG1 

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022.
% See file license.txt in installation directory for licensing terms.

function c = asymp_dtf(u,A,pf,nFreqs,metric,alpha)

if ~(nargin == 6)
   error('ASYMP_DTF requires six input arguments.')
end
[m,n] = size(u);
if m > n
   u = u.';
end
np = length(u);
[nChannels,~,p] = size(A);
Af = A_to_f(A, nFreqs);

flgVerbose = 0;

% Variables pre-alocation
dtf  = zeros(nChannels,nChannels,nFreqs); % dtf will hold complex DTF estimates
dtf2 = zeros(nChannels,nChannels,nFreqs);

if alpha ~= 0
   th  = zeros(nChannels,nChannels,nFreqs);
   ci1 = zeros(nChannels,nChannels,nFreqs);
   ci2 = zeros(nChannels,nChannels,nFreqs);
   varass1 = zeros(nChannels,nChannels,nFreqs);
   varass2 = zeros(nChannels,nChannels,nFreqs);
   patdfr = zeros(nChannels,nChannels,nFreqs);
   patdenr = zeros(nChannels,nChannels,nFreqs);
   pvalues = zeros(nChannels,nChannels,nFreqs);

   if flgVerbose
      switch lower(metric)
         case {'euc'}
            disp('* Original DTF and asymptotic statistics')
         case {'diag'}
            disp('* Generalized DTF or DC and asymptotic statistics')
         case {'info'}
            disp('* Information DTF and asymptotic statistics')
         otherwise
            error('Unknown metric.')
      end
   end
elseif flgVerbose
   switch lower(metric)
      case {'euc'}
         disp('* Original DTF estimation')
      case {'diag'}
         disp('* Generalized DTF/DC estimation')
      case {'info'}
         disp('* Information DTF estimation')
      otherwise
         error('Unknown metric.')
   end
end

switch lower(metric)
   case {'euc'}               % for DTF
      ddtf_dev = zeros(1,nChannels^2);
      pfe  = eye(nChannels);  % for complex DTF calculation

   case {'diag'}              % for DC
      evar_d = mdiag(pf);
      evar_d_big = kron(eye(2),kron(evar_d,eye(nChannels)));

      debig_de = fdebig_de_dtf(nChannels); %New Theta_K
      dedinv_deh = debig_de*diag(vec(eye(nChannels)));
      pfe = pf;                % for complex DC calculation

   case {'info'}               % for iDTF
      evar_d = mdiag(pf);
      evar_d_big = kron(eye(2),kron(evar_d,eye(nChannels)));
      evar_big = kron(eye(2),kron(pf,eye(nChannels)));

      debig_de = fdebig_de_dtf(nChannels); %New Theta_K
      dedinv_deh = debig_de*diag(vec(eye(nChannels)));
      pfe = pf;                % for complex iDTF calculation

   otherwise
      error('Unknown metric.')
end

gamma = bigautocorr(u, p);
omega = kron(pinv(gamma), pf);
omega_evar = 2*Dup(nChannels)*pinv(Dup(nChannels))*kron(pf, pf) ...
                            *(pinv(Dup(nChannels)).')*Dup(nChannels).';

icdf_norm_alpha = norminv(1 - alpha/2.0,0,1);

for ff = 1:nFreqs
   f = (ff - 1)/(2*nFreqs); % Corrected 7/25/2011, f starts at 0.
   Ca = fCa(f, p, nChannels);

   Af_ff = reshape(Af(ff,:,:),[nChannels, nChannels]);
   Hf = pinv(Af_ff); h = Hf(:);  % Equivalent to h = vec(Af[ff, :, :].I)

   h = [real(h); imag(h)];    % h = cat(h.real, h.imag, 0)
   H = fdh_da(Af_ff);         % = ha; H = fdh_da(mat(Af[ff, :, :]), n)

   Omega_h = H*Ca*omega*Ca.'*H.';  % \Omega_h
   L = fChol(Omega_h); % real-part only

   for i = 1:nChannels
      Ii = fIi(i,nChannels);

      switch lower(metric)
         case {'euc'}           % for DTF
            Iie  = Ii;
         case {'diag'}          % for DC
            Iie  = Ii*evar_d_big*Ii;
         case {'info'}          % for iDTF
            Iie  = Ii*evar_big;
         otherwise
            error('Unknown metric.')
      end

      for j = 1:nChannels
         Iij = fIij(i,j,nChannels);

         switch lower(metric)
            case {'euc'}               % for DTF
               Iije = Iij;

            case {'diag'}              % for DC
               Iije = Iij*evar_d_big*Iij;

            case {'info'}              % for iDTF
               Iije = Iij*evar_d_big;

            otherwise
               error('Unknown metric.')
         end

         num = h.'*Iije*h;
         den = h.'*Iie*h;

         dtf2(i,j,ff) = num/den; % |DTF_{ij}(ff)|^2 squared-|DTF|
         dtf(i,j,ff)  = Hf(i,j)*sqrt(pfe(j,j))/sqrt(den); % complex-DTF

         if alpha ~= 0
            %'Add evar differentiation'
            switch lower(metric)
               case {'euc'}               % for DTF
                  %nop

               case {'diag'}              % for DC
                  dnum_dev = kron((Iij*h).', h.'*Iij)*dedinv_deh;
                  %'derivative of den by vecE'
                  dden_dev = kron((Ii*h).', h.'*Ii)*dedinv_deh;
                  ddtf_dev = (den*dnum_dev - num*dden_dev)/(den^2);

               case {'info'}              % for iDTF
                  %'derivative of num by vecE'
                  dnum_dev = kron((Iij*h).', h.'*Iij) * dedinv_deh;
                  %'derivative of den by vecE'
                  dden_dev = kron((Ii*h).', h.'*Ii) * debig_de;
                  ddtf_dev = (den*dnum_dev - num*dden_dev)/(den^2);

               otherwise
                  error('Unknown metric.')
            end

            G1h = 2*h.'*Iije/den - 2*num*h.'*Iie/(den^2); % Eq. (15)
            %   G1 = -G1h*H*Ca;                           % (Cont.)
            varalpha = G1h*Omega_h*G1h.';
            varevar = ddtf_dev*omega_evar*ddtf_dev.';
            varass1(i,j,ff) = (varalpha + varevar)/np;    % Eq. (14)

            ci1(i,j,ff) = dtf2(i,j,ff) ...
                          - sqrt(varass1(i,j,ff))*icdf_norm_alpha;
            ci2(i,j,ff) = dtf2(i,j,ff) ...
                          + sqrt(varass1(i,j,ff))*icdf_norm_alpha;

            G2h = Iije/den;
            d = fEig(real(L), real(G2h)); % real() 28May2013

            patdf = (sum(d).^2)./sum(d.^2);
            patden = sum(d)./sum(d.^2);

             th(i,j,ff) = chi2inv((1 - alpha), patdf)/(patden*np); % original KS
%            th(i,j,ff) = chi2inv((1 - alpha), patdf)./(np);
            pvalues(i,j,ff) = 1 - chi2cdf(dtf2(i,j,ff)*patden*np, patdf);

            varass2(i,j,ff) = patdf/(patden*np).^2;
            patdfr(i,j,ff) = patdf;
            patdenr(i,j,ff) = patden;

         else % as alpha == 0, do not compute asymptotics
            %nop
         end
      end
   end
end

if alpha ~= 0
   c.dtf = dtf;   % Complex DTF2/DC2/iDTF2
   c.dtf2 = dtf2;
   c.th = th;
   c.ci1 = ci1;  % Lower CI
   c.ci2 = ci2;  % Upper CI
   c.metric = metric;
   c.alpha = alpha;
   c.p = p;
   c.pvalues = pvalues; % p-values associated to DTF2/DC2/iDTF2
   c.patden = patdenr;
   c.patdf = patdfr;
   c.varass1 = varass1;
   c.varass2 = varass2;

   % Statistically significant DTF2 on frequency scale
   dtf2_temp = ((abs(dtf2) - abs(th)) > 0).*dtf2 ...
                                           + ((abs(dtf2) - abs(th)) <= 0)*(-1);
   dtf2_temp(dtf2_temp < 0) = NaN; % Octave
   c.dtf2_th = dtf2_temp;

else
    c.dtf = dtf;
    c.dtf2 = dtf2;
    c.metric = metric;
    c.alpha = 0;
    c.p = p;
    c.th = [];
    c.ci1 = [];
    c.ci2 = [];
    c.pvalues = [];
    c.patden = [];
    c.patdf = [];
    c.varass1 = [];
    c.varass2 = [];
end

% Power spectra and spectral coherence calculation
c.SS = ss_alg(A, pf, nFreqs);
c.coh2 = coh_alg(c.SS);
end

%==========================================================================
function gamma = bigautocorr(x, p)
% Autocorrelation. Data in rows. From order 0 to p-1.
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
% Returns the eigenvalues 

%L = mat(cholesky(omega, lower=1))
D = L.'*G2*L;
%    d = eigh(D, eigvals_only=True)
%disp('fEig: eig or svd?')
d = svd(D);
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
% Returns Iij of the formula
Iij = zeros(1,n^2);
Iij(n*(j - 1) + i) = 1;
Iij = diag(Iij);
c = kron(eye(2), Iij);
end

%==========================================================================
function c = fIi(i,n)
% Returns Ii of the formula 
Ii = zeros(1,n);
Ii(i) = 1;
Ii = diag(Ii);
Ii = kron(eye(n), Ii);
c = kron(eye(2), Ii);
end

%==========================================================================
function d = fCa(f, p, n)
% Returns C* of the formula 
C1 = cos(-2*pi*f*(1:p));
S1 = sin(-2*pi*f*(1:p));
C2 = [C1; S1];
d = kron(C2, eye(n^2));
end

%==========================================================================
function t = TT(a,b)
%  TT(a,b)*vec(B) = vec(B.T), where B is (a x b). 
t = zeros(a*b);
for i = 1:a
   for j =1:b
      t((i - 1)*b + j,(j - 1)*a + i) = 1;
   end
end
t = sparse(t);
end

%==========================================================================
function c = fdebig_de_dtf(n)
%  New \Theta_K for DTF asymptotics 
%c = kron(kron(eye(2*n),TT(n, 2*n)), eye(n)));
%A = sparse(kron(eye(2*n), TT(n, 2*n)));
%c = sparse(kron(A, eye(n)));
A  = sparse(kron(TT(n^2,2),eye(n^2)));
A1 = sparse(kron(eye(2),A));

% A=sparse(kron(TT(n^2,1),eye(n)));  % The code line corrected by
A  = sparse(kron(TT(n,n),eye(n)));   % Rezaei et. al. (2022).

A2 = sparse(kron(eye(n),A));
A3 = sparse(kron(eye(n^2),vec(eye(n))));
A4 = sparse(A2*A3);
A5 = sparse(kron(vec(eye(2)),A4));
c  = sparse(A1*A5);
end

%==========================================================================
function c = vec(x)
% vec = lambda x: mat(x.ravel('F')).T
c = x(:);
end

%==========================================================================
function L = fChol(omega)
% Try Cholesky factorization
try
   L = chol(omega)';
   % If there's a small negative eigenvalue, diagonalize
catch % err
   %   disp('linalgerror, probably IP = 1.')
   [v,d] = eig(omega);
   L = zeros(size(v));
   for i =1:length(d)
      if d(i,i)<0
         d(i,i)=eps;
      end
      L(:,i) = v(:,i)*sqrt(d(i,i));
   end
end
end

%==========================================================================
function c = mdiag(a)
%  diagonal matrix
c = diag(diag(a));
end

%==========================================================================
function d = Dup(n)
% D*vech(A) = vec(A), with symmetric A 
d = zeros(n*n, (n*(n + 1))/2);
count = 1;
for j= 1:n
   for i = 1:n
      if i >= j
         d((j - 1)*n + i,count)=1;
         count = count + 1;
      else
         d((j - 1)*n + i,:)=d((i - 1)*n + j,:);
      end
   end
end
end

%==========================================================================
function hh=fdh_da(Af)
% Derivative of vec(H) by vec(A), with $H = A^{-1}$ and complex A. 
ha = pinv(Af);
h  = -kron(ha.', ha);

h1 = [real(h) -imag(h)];
h2 = [imag(h) real(h)];   %h2 = cat(h.imag, h.real, 1)
hh = -[h1; h2];           %hh = cat(h1, h2, 0)
end
%==========================================================================

%%
%        1         2         3         4         5         6         7         8         9
%23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
%% Change Log:
% [2011/07/25]: The asymp_pdc routine, which asymp_dtf is derived from, was corrected on
%              to match the frequency range with plotting routine, f = 0 was
%              included in the "frequency" for-loop:
%                                for ff = 1:nFreqs,
%                                   f = (ff-1)/(2*nFreqs); %
%                                        ^?^^
% [2021/09/02, LAB]: complex dtf calculation added
%
% [2015/01/07]: Optimization \(^o^)/
%
% [2022/11/05]: Line 339: Correction in the function c = fdebig_de_dtf(n)
%               pointed out by Rezaei et al. (2022). Fast Asymptotic
%               Algorithm for Real-Time Causal Connectivity Analysis of
%               Multivariate Systems and Signals. Signal Processing.
%               https://doi.org/10.1016/j.sigpro.2022.108822

% [EOF]
