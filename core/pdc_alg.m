%% PDC_ALG
%        Compute partial directed coherence measure given by metric "option"
%        from series j--> i.
%
%% Syntax:
%        c = PDC_ALG(u,alg,criterion,nFreqs,metric,maxIP,alpha)
%
%% Input arguments:
%        u         - data
%        alg       - algorithm (1: Nutall-Strand),(2: mlsm) ,
%                              (3: Vieira Morf),  (4: QR artfit)
%        criterion - AR order selection criteria =>
%                                   1: AIC; 2: Hanna-Quinn; 3: Schwarz;
%                                   4: FPE, 5: fixed order in MaxIP
%        nFreqs - number of point in [0,fs/2] frequency scale
%        metric   'euc'  - Euclidean ==> original PDC
%                 'diag' - diagonal ==> gPDC (generalized )
%                 'info' - information ==> iPDC
%        maxIP - externally defined maximum IP
%        alpha - Significance level for null hypothesis testing
%                alpha = .05 is default;
%                if alpha = zero, no asymptotic statistics is computed.
%
%% Output arguments:
%        c structure variable with following fields:
%        |-- .pdc       - complex PDC estimates
%        |-- .pdc2      - |PDC(f)|^2 estimate
%        |-- .pdc2_th   - above threshold |PDC(f)|^2 otherwise PDC2(f) = NaN.
%        |-- .pvalues   - p-values associated to PDC2 estimates. 
%        |-- .th        - |PDC|^2 threshold value with (1-avalue) significance level.
%        |-- .{ic1,ic2} - upper and lower (1 - alpha) confidence interval of |PDC|^2 estimates
%        |-- .metric    - metric used for PDC calculation 
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
%        c.{pdc,pdc2,pvalues,th,pdc2_th,ic1,ic2,metric,alpha,p,patdenr,patdfr,SS,coh2}
%
%% Example:
%                 u=sunmeladat([4 3]);  % Andrews & Herzberg 1936-1972
%                                       % sunspot-melanoma series
%                 u=detrend(u);         % Detrend the series 
%                 c=pdc_alg(u,64,'diag',1,1,30,0.01); 
%                 figure; xplot(c); %pretty plot
%
%% See also: ASYMP_PDC, DTF_ALG

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.

function c = pdc_alg(u,nFreqs,metric,alg,criterion,maxIP,alpha)

% If the number of input parameters is smaller than seven, following default
% values are assumed for PDC calculation.
if nargin < 7, alpha = 0; end;      % do not calculate asymptotic statistics
if nargin < 6, maxIP = 30; end;     % defaults value 
if nargin < 5, criterion =  1; end; % AIC order choice
if nargin < 4, alg = 1;  end;       % Nutall-Strand is default AR estimator
if nargin < 3, metric = 'diag'; end;% gPDC estimation
if nargin < 2, nFreqs = 128; end;   % 128 points on freq scale.
 
[m,n] = size(u);
if m > n,
   u = u.';
end;

nSegLength = length(u);

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);

%==========================================================================
%    Testing for adequacy of MAR model fitting through Portmanteau test
%==========================================================================
h = 20; % testing lag
VARadequacy_signif = 0.05;
aValueVAR = 1 - VARadequacy_signif;
flgPrintResults = 1;
[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueVAR,h,...
                                                          flgPrintResults);
%==========================================================================
%            PDC, threshold and confidence interval calculation.
%==========================================================================

% if alpha == 0, no asymptotic statistics calculation is performed, and ASYMP_PDC
% returns only the PDC. This option is much faster!!
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha);

c.A = A; 
c.pf = pf; 
c.nfreqs = nFreqs;
c.criterion = criterion; 
c.alg = alg;
c.Pass = Pass; 
c.Portmanteau = Portmanteau; 

% Power spectra and coherence calculation
c.SS = ss_alg(A, pf, nFreqs);
c.coh = coh_alg(c.SS);

% Statistically significant PDC on frequency scale
if alpha ~= 0
   pdc2_temp = ((abs(c.pdc2)-c.th) > 0).*c.pdc2 + ((abs(c.pdc2)-c.th) <= 0)*(-1);
   pdc2_temp(ind2sub(size(pdc2_temp),find(pdc2_temp == -1))) = NaN;
   c.pdc2_th = pdc2_temp;
else
   c.pdc2_th = [];
end
end