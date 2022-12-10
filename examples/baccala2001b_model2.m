%% BACCALA & SAMESHIMA (2001b) - Model 2: 6-dim VAR[4] Scalp Model
%
% Description:
%
% Baccala & Sameshima. Overcoming the limitations of correlation analysis 
% for many simultaneously processed neural structures, Progress in Brain 
% Research, 130:33--47, 2001.
%
% <http://dx.doi.org/10.1016/S0079-6123(01)30004-3>
% 
% Model II: 6-dimension VAR(4) Scalp model
%
%% See also: mvar, mvarresidue, asymp_pdc, asymp_dtf, gct_alg, 
%              igct_alg, xplot, xplot_pvalues             

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.

clear; clc; format compact; format short

%% Data sample generation

% Seeding random number generator
rng('default')
%rng('shuffle')

nBurnIn = 5000;    % number of points discarded at beginning of simulation
nPoints  = 2000;   % number of analyzed samples points
N=nBurnIn+nPoints; % number of simulated points

u = fbaccala2001b_model2( nPoints, nBurnIn );

[nSegLength,nChannels]=size(u);
chLabels = {'x_1';'x_2';'x_3';'x_4';'x_5';'x_6'}; % or chLabels = [];
fs = 1;

%% Interaction diagram
%
% <<fig_baccala2001b_model2_graph.png>>
%
% Figure 9 from Baccala & Sameshima. _Progr in Brain Res_, 130:33--47, 2001.

%% Equation
%
% <<fig_baccala2001b_model2_eq.png>>
%

%%
% Data checking and pre-processing: detrending and normalization options

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

disp('Running MVAR estimation routine.')
maxIP = 30;         % maximum model order to consider.
alg = 1;            % 1: Nutall-Strand MVAR estimation algorithm
criterion = 1;      % 1: AIC, Akaike Information Criteria
[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);

disp(['Number of channels = ' int2str(nChannels) ' with ' ...
    int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%%
% Testing for adequacy of MAR model fitting through Portmanteau test
%

h = 20; % testing lag
MVARadequacy_signif = 0.05; % VAR model estimation adequacy significance
                            % level
aValueMVAR = 1 - MVARadequacy_signif; % Confidence value for the testing

flgPrintResults = 1; % Print results on Command Window
[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueMVAR,h,...
                                           flgPrintResults);

%% Granger causality test (GCT) and instantaneous GCT
%

gct_signif  = 0.01;  % Granger causality test significance level
igct_signif = 0.01;  % Instantaneous GCT significance level
flgPrintResults = 1;

[Tr_gct, pValue_gct]   =  gct_alg(u,A,pf, gct_signif,flgPrintResults);
[Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgPrintResults);


%% Original PDC estimation
%
% PDC analysis results are saved in *c* struct.
% See asymp_pdc.m for more detail. 

nFreqs = 128;  % Number of points in the frequency scale.
alpha  = 0.01; % Significance level for PDC testing 
metric = 'euc';  % euc  = Euclidean, original PDC or DTF;
                 % diag = generalized PDC, gPDC, or  gDTF / DC;
                 % info = information PDC, iPDC, or  iDTF.
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
c.Tragct = Tr_gct;
c.pvaluesgct = pValue_gct;

%% $|PDC(\lambda)|^2$ Matrix-Layout Plotting

flgPrinting = [1 1 1 2 3 1 2];
flgColor = 1;
flgScale = 1;
flgMax = 'tci';
flgSignifColor = 3; % Line color for PDC2/DTF2 plotting on the frequency scale.
%                     3: red (significat) / green (not significant) 
w_max = fs/2;

strTitle = ['6-dimension linear VAR[4] Model II: [N=' int2str(nSegLength) ...
            'pts; IP=' int2str(c.p) ']'];
strID = 'Baccala & Sameshima (2001b) Model II'; % figure window bar ID
[h1,~,~] = xplot(strID,c,flgPrinting,fs,w_max,chLabels, ...
                                       flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric,'pdc', strTitle);

%% 
% Result from the original article, Baccala & Sameshima (2001b) 
% Figure 6b.
%
% <<fig_baccala2001b_model2_pdc_result.png>>
% 

%% Original DTF estimation
% DTF analysis results are saved in *d* struct.
%

d = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Estimate DTF and asymptotic statistics

%% $|DTF(\lambda)|^2$ Matrix-Layout Plotting
%

flgPrinting = [1 0 1 2 0 1 2]; % not plotting threshold and coherence
[h2,~,~] = xplot(strID,d,flgPrinting,fs,w_max,chLabels, ...
                                       flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric,'dtf', strTitle);


%% Animation
% Gif
% <<fig_gif_baccala2001b_escalp_eeg_pdc.gif>>


