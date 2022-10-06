%% EICHLER ET AL.(2006) - VAP II 5-dimension VAR[3] process
%
% An example borrowed from Schelter et al. (2009) 
%
% Schelter, Timmer, Eichler. Assessing the strength of directed influences
% among neural signals using renormalized PDC. J Neurosci Methods,
% 179:121-130, 2009.
%   [http://dx.doi.org/10.1016/j.jneumeth.2009.01.006]
% 
% 3.2 Vector autoregressive process II (Eqs. 16-20, page 124)
%
%% See also: mvar, mvarresidue, asymp_pdc, asymp_dtf, gct_alg, 
%              igct_alg, xplot, xplot_pvalues             

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.

clear; clc

%% Equation Model II with feedback
% 3.2 Vector autoregressive process II (Eqs. 16-20, page 124)
%
% <<fig_schelter2009_vap2_equations.png>>
%


%% Data generation
%

nDiscard = 10000;    % number of points discarded at beginning of simulation
nPoints  = 3000;   % number of analyzed samples points

u = fschelter2009_vap2(nPoints, nDiscard);
%chLabels = []; % or 
chLabels = {'x_1';'x_2';'x_3';'x_4';'x_5'};
fs =1;

%%
% Data pre-processing: detrending and normalization options
flgDetrend = 1;     % Detrending the data set
flgStandardize = 0; % No standardization

[nChannels,nSegLength]=size(u);
if nChannels > nSegLength, u=u.'; 
   [nChannels,nSegLength]=size(u);
end

if flgDetrend
   for i=1:nChannels, u(i,:)=detrend(u(i,:)); end
   disp('Time series were detrended.');
end

if flgStandardize
   for i=1:nChannels, u(i,:)=u(i,:)/std(u(i,:)); end
   disp('Time series were scale-standardized.');
end

%% MVAR model estimation

maxIP = 30;         % maximum model order to consider.
alg = 1;            % 1: Nutall-Strand MVAR estimation algorithm;
%                   % 2: minimum least squares methods;
%                   % 3: Vieira Morf algorithm;
%                   % 4: QR ARfit algorith.

criterion = 1;      % Criterion for order choice:
%                   % 1: AIC, Akaike Information Criteria; 
%                   % 2: Hanna-Quinn;
%                   % 3: Schwarz;
%                   % 4: FPE;
%                   % 5: fixed order given by maxIP value.

disp('Running MVAR estimation routine...')

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);

disp(['Number of channels = ' int2str(nChannels) ' with ' ...
    int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%==========================================================================
%    Testing for adequacy of MAR model fitting through Portmanteau test
%==========================================================================
h = 20; % testing lag
MVARadequacy_signif = 0.05; % VAR model estimation adequacy significance
                            % level
aValueMVAR = 1 - MVARadequacy_signif;
flgPrintResults = 1;
[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueMVAR,h,...
                                           flgPrintResults);

%% Granger causality test (GCT) and instantaneous GCT

gct_signif  = 0.05;  % Granger causality test significance level
igct_signif = 0.05;  % Instantaneous GCT significance level
flgPrintResults = 1; % Flag to control printing gct_alg.m results on command window.
[Tr_gct, pValue_gct] = gct_alg(u,A,pf,gct_signif,flgPrintResults);
[Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgPrintResults);

%% Information PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_pdc.m or issue 
%
%   >> help asymp_pdc 
%
% command for more detail.

nFreqs = 128;
metric = 'info';  % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
alpha = 0.05;

c = asymp_pdc(u,A,pf,nFreqs,metric,alpha);
c.Tragct = Tr_gct;         % Assigning GCT results to c struct variable.
c.pvaluesgct = pValue_gct;

%% iPDC2 Matrix-Layout Plot

flgPrinting = [1 1 1 2 3 0 2]; % overriding default setting
flgColor = 1;
w_max=fs/2;
alphastr = int2str(100*alpha);

strID = 'Schelter et al. (2009)';
strTitle = ['Schelter et al. (2009) linear model II: ' ...
                               int2str(nPoints) ' data points.'];
[h,~,~] = xplot(strID,c,flgPrinting,fs,w_max,chLabels,flgColor,2,'pdc');
xplot_title(alpha,metric,'pdc',strTitle);


%% Result from Schelter et al.(2009) 
% Figure 2, page 125.
%
% <<fig_schelter2009_vap2_result.png>>
% 

%% 
% Figure caption: 
% "Fig 2. Renormalized partial directed coherence (off-diagonal) for the
% example of the VAR[2] process. The results are sorted as a matrix, where
% in the ith row and the jth column the influence from process j onto
% process i is displayed. The black line represents the renormalized
% partial directed coherence values while the gray areas mark the
% corresponding 95% confidence intervals of one single realization. From
% 1000 realizations empirical confidence intervals have been derived. These
% are indicated by the dashed black lines. Only the valid interactions are
% revealed by this analysis." (Reproduced from  Schelter et al. J Neurosci
% Methods, 179:121-130, 2009)

%% Remarks:
%
% * Our $|_iPDC(\lambda)|^2$ results are from a single simulation data
% while normalized PDC are from 1000 simulations. Note that the y-axis
% scales for normalized PDC by Schelter et al. (2009) are different from
% iPDC2 plots.
%
% * Note again that for linear model the mean amplitude of PDC estimates
% is roughly proportional to  relative coefficient values of ') the
% autoregressive model.
%
% * One may see occasional false-positive inferences in iPDC estimates for
% $\alpha = 5\%$  null-hypothesis testing, however as one might notice the
% magnitude should most likely (from our experience) be smaller than 0.01
% and Granger causality test p-values printed just above each subplot
% should be close to chosen $\alpha$ values.
