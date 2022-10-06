%% BACCALA & SAMESHIMA (2001b) - Model 2: 6-dim VAR[4] Scalp Model variant
%
% Description:
%
% Variant of model borrowed from 
% Baccala & Sameshima. Overcoming the limitations of correlation analysis 
% for many simultaneously processed neural structures, Progress in Brain 
% Research, 130:33--47, 2001.
%
% <http://dx.doi.org/10.1016/S0079-6123(01)30004-3>
% 
% 6-dimensional VAR(4) Scalp variant Model II 
%
%% See also: mvar, mvarresidue, asymp_pdc, asymp_dtf, gct_alg, 
%              igct_alg, xplot, xplot_pvalues             

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.

clear; clc; format compact; format short

%% Interaction diagram
%
% <<fig_baccala2001b_model2_variant_graph.png>>
%
% Figure 2a from Baccala & Sameshima. _Biol. Cybern._ *84*:463-474, 2001.

%% Equation
%
% <<fig_baccala2001b_model2_variant_eq.png>>
%

%% Data sample generation

nDiscard = 5000;    % number of points discarded at beginning of simulation
nPoints  = 1000;    % number of analyzed samples points
N=nDiscard+nPoints; % number of simulated points

u = fbaccala2001b_model2_variant( nPoints, nDiscard );

chLabels = {'x_1';'x_2';'x_3';'x_4';'x_5';'x_6'};

fs = 1;

%% Data pre-processing: detrending and normalization options
%

flgDetrend = 1;     % Detrending the data set
flgStandardize = 0; % No standardization
[nChannels,nSegLength] =size(u);
if nChannels > nSegLength
   u = u.'; 
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
%

maxIP = 30;         % maximum model order to consider.
alg = 1;            % 1: Nutall-Strand MVAR estimation algorithm
criterion = 1;      % 1: AIC, Akaike Information Criteria

disp('Running MVAR estimation routine.')

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);

disp(['Number of channels = ' int2str(nChannels) ' with ' ...
    int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%%
% Testing for adequacy of MAR model fitting through Portmanteau test

h = 20; % testing lag
MVARadequacy_signif = 0.05; % VAR model estimation adequacy significance
                            % level
aValueMVAR = 1 - MVARadequacy_signif; % Confidence value for the testing
flgPrintResults = 1;

[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueMVAR,h,...
                                           flgPrintResults);

%% Granger causality test (GCT) and instantaneous GCT
%

gct_signif  = 0.001;  % Granger causality test significance level
igct_signif = 0.001;  % Instantaneous GCT significance level

flgPrintResults = 1;

[Tr_gct, pValue_gct]   =  gct_alg(u,A,pf, gct_signif,flgPrintResults);
[Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgPrintResults);


%% Original PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_pdc.m.

metric = 'euc';  % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
nFreqs = 128;
alpha = 0.001;
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
c.pvaluesgct = pValue_gct; % Necessary for printing GCT
c.Tragct = Tr_gct;

%% $|PDC(\lambda)|^2$ Matrix Layout Plotting
%

flgPrinting = [1 1 1 2 3 1 2]; % overriding default setting
flgColor = 0;
w_max=fs/2;

strTitle1 = ['6-dimensional linear VAR[4] Variant Model II: '];
strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p) ']'];
strTitle =[strTitle1 strTitle2];

strWindowBar = 'Baccala & Sameshima (2001b) Model II Variant';
[h,~,~] = xplot(strWindowBar,c,...
                          flgPrinting,fs,w_max,chLabels,flgColor,2,'all');
xplot_title(alpha,metric,'pdc', strTitle);

