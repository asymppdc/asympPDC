%% DTF_ALG
%        Compute directed transfer function measure given by the metric option
%        from series j--> i.
%
%% Syntax:
%        c = dtf_alg(u,alg,criterion,nFreqs,metric,maxIP,alpha)
%
%% Input arguments:
%        u         - data
%        alg       - algorithm (1: Nutall-Strand); (2: mlsm);
%                              (3: Vieira Morf); (4: QR artfit)
%        criterion - AR order selection criteria:
%                                   1: AIC; 2: Hanna-Quinn; 3: Schwarz;
%                                   4: FPE, 5: fixed order in MaxIP
%        nFreqs - number of point in [0,fs/2] frequency scale
%        metric - 'euc'  - Euclidean ==> DTF;
%                 'diag' - diagonal ==> DC;
%                 'info' - information ==> iDTF.
%        maxIP - externally defined maximum IP
%        alpha - Significance level for null hypothesis testing,
%                alpha = .05 is default; if alpha = zero, asymptotic statistics
%                is not computed.
%
%% Output arguments:
%        c struct variable with following fields:
%        |-- .dtf       - complex DTF estimates
%        |-- .dtf2      - |DTF|^2 estimates
%        |-- .dtf2_th   - Statistically significant |DTF|^2 estimates on frequency scale
%        |-- .pvalues   - p-values associated to DTF2 estimates. 
%        |-- .th        - |DTF|^2 threshold value with (1-alpha) significance level.
%        |-- .{ic1,ic2} - upper and lower (1 - alpha) confidence interval of |DTF|^2 estimates
%        |-- .metric    - metric used for DTF calculation 
%        |-- .alpha     - significance level
%        |-- .criterion - AR order selection criterion
%        |-- .alg       - AR estimation algorithm
%        |-- .A         - AR coefficient matrix estimate
%        |-- .pf        - covariance matrix estimate
%        |-- .p         - VAR model order
%        |-- .patdenr   - 
%        |-- .patdfr    - degree of freedom 
%        |-- .SS        - power spectra
%        +-- .coh2      - squared spectral coherence
%    or
%        c.{dtf,dtf2,pvalues,th,ic1,ic2,metric,alpha,p,patdenr,patdfr,SS,coh2}
%
%% Example:
%                 u = sunmeladat([4 3]);  % Andrews & Herzberg 1936-1972
%                                         % sunspot-melanoma series
%                 u = detrend(u);         % Detrend the series 
%                 c = dtf_alg(u,64,'diag',1,1,30,0.01); 
%                 figure; xplot(c); %pretty plotting     
%
%% See also: ASYMP_DTF, PDC_ALG

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.

%%

function c = dtf_alg(u,nFreqs,metric,alg,criterion,maxIP,alpha)

% If the number of input parameters is smaller than seven, following default
% values are adopted for DTF calculation.
if nargin < 7, alpha = 0; end       % do not calculate asymptotic statistics
if nargin < 6, maxIP = 30; end      % defaults value 
if nargin < 5, criterion =  1; end  % AIC order choice
if nargin < 4, alg = 1;  end        % Nuttall-Strand is default AR estimator
if nargin < 3, metric = 'diag'; end % DC estimation
if nargin < 2, nFreqs = 128; end    % 128 points on the frequency scale.
 
[m,n] = size(u);
if m > n
   disp('Data vector seems to be transposed.')
   u = u.';
end

nSegLength = length(u);

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);

%==========================================================================
%    Testing the adequacy of MAR model fitting through Portmanteau test
%==========================================================================
h = 20; % testing lag
VARadequacy_signif = 0.05;
aValueVAR = 1 - VARadequacy_signif;
flgPrintResults = 1;
[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueVAR,h,...
                                                          flgPrintResults);
%==========================================================================
%            DTF, threshold and confidence interval calculation.
%==========================================================================

% if alpha == 0, no asymptotic statistics is performed. ASYMP_DTF returns
% only the DTF. This option is much faster!!
c = asymp_dtf(u,A,pf,nFreqs,metric,alpha);

c.A = A; 
c.pf = pf; 
c.nfreqs = nFreqs;
c.criterion = criterion; 
c.alg = alg;
c.Pass = Pass; 
c.Portmanteau = Portmanteau; 

% Power spectra and spectral coherence calculation
c.SS = ss_alg(A, pf, nFreqs);
c.coh2 = coh_alg(c.SS);

% Statistically significant DTF on frequency scale
if alpha ~= 0
   dtf2_temp = ((abs(c.dtf2)-c.th) > 0).*c.dtf2 + ((abs(c.dtf2)-c.th) <= 0)*(-1);
   dtf2_temp(ind2sub(size(dtf2_temp),find(dtf2_temp == -1))) = NaN;
   c.dtf2_th = dtf2_temp;
else
   c.dtf2_th = [];
end
