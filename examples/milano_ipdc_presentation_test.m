%% Baccala & Sameshima (2001b) 7-dimension VAR[2] model with loop and feedback
%
%% Description:
%
% Baccala & Sameshima. Overcoming the limitations of correlation analysis 
% for many simultaneously processed neural structures, Progress in Brain 
% Research, 130:33--47, 2001.
%
% <https://dx.doi.org/10.1016/S0079-6123(01)30004-3>
%
%% See also:
% 
% # Koichi Sameshima, Daniel Y. Takahashi, Luiz A. Baccala. On the
% statistical performance of Granger-causal connectivity estimators. Brain
% Informatics (2015) 2:119?133.
%
% <https://dx.doi.org/10.1007/s40708-015-0015-1>
% 
% Example Model 1 - 7-dimension VAR[2] model with loop and feedback
%
%% Other routines
%  See also  fbaccala2001b_model1_feedback
% < milano_ipdc_presentation_test.html |milano_ipdc_presentation_test|>

% (C) Koichi Sameshima & Luiz A. Baccala', 2022. See file license.txt in
% installation directory for licensing terms.


%% Data sample generation

clear; clc; format compact; format short

nDiscard = 50000;          % discarded points at beginning of simulation
nSegLength  = 5000;        % number of analyzed samples points
N = nDiscard + nSegLength; % number of simulated points

u = fbaccala2001b_model1_feedback( nSegLength, nDiscard );

[nSegLength,nChannels] = size(u);

%chLabels = []; % or 
chLabels = {'x_1';'x_2';'x_3';'x_4';'x_5';'x_6';'x_7'};

fs = 1;
%% Model Interaction diagram, model equation and expected results
%
% <<fig_baccala2001b_graph.png>>
%
% Figure 2a from Baccala & Sameshima. _Biol. Cybern._ *84*:463-474, 2001.

%% Equation Model I with feedback
%
% <<fig_baccala2001b_eq.png>>
%

%% iDTF and iPDC expected results
%
% <<fig_sameshima2005_fig2_iDTFiPDCresults.png>>
%

%% Data pre-processing: detrending and standardization options
%

flgDetrend = 1;     % Detrending the data set
flgStandardize = 0; % No standardization
[nChannels,nSegLength] =size(u);
if nChannels > nSegLength, 
   u = u.'; 
   [nChannels,nSegLength]=size(u);
end;
if flgDetrend,
   for i=1:nChannels, u(i,:)=detrend(u(i,:)); end;
   disp('Time series were detrended.');
end;
if flgStandardize,
   for i=1:nChannels, u(i,:)=u(i,:)/std(u(i,:)); end;
   disp('Time series were scale-standardized.');
end;

%% 
% MVAR model estimation

alpha = 0.01;

maxIP = 30;         % maximum model order to consider.
alg = 1;            % 1: Nutall-Strand MVAR estimation algorithm
criterion = 1;      % 1: AIC, Akaike Information Criteria

fprintf('\n')
disp('Running MVAR estimation routine.')

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);

disp(['Number of channels = ' int2str(nChannels) ' with ' ...
    int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);
   fprintf('\n')

%%
% Testing for adequacy of MAR model fitting through Portmanteau test

h = 20; % testing lag
MVARadequacy_signif = 0.05; % VAR model estimation adequacy significance
                            % level
aValueMVAR = 1 - MVARadequacy_signif; % Confidence value for the testing

flgPrintResults = 1;

[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueMVAR,h,...
                                           flgPrintResults);

%%
% Granger causality test (GCT) and instantaneous GCT (iGCT)

gct_signif  = alpha;  % Granger causality test significance level
igct_signif = alpha;  % Instantaneous GCT significance level

metric = 'info'; % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
flgPrintResults = 1;

[Tr_gct, pValue_gct] = gct_alg2(u,A,pf,gct_signif,flgPrintResults);
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

c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
    c.Tragct = Tr_gct;
    c.pvaluesgct = pValue_gct;

%%
%
% $|_iPDC(\lambda)|^2$ Matrix Layout Plotting
%

flgPrinting = [1 0 1 2 3 0 2]; % overriding default setting
flgScale = 1;
flgMax = 'pdc';
flgSignifColor = 3;
flgColor = [0];

w_max=fs/2;

strTitle1 = ['Model I with feedback: '];
strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p)  ' ]'];
strTitle =[strTitle1 strTitle2];
strID = 'Baccala & Sameshima (2001) Example 3';

[hxlabel,hylabel] = xplot(strID,c,flgPrinting,fs,w_max,chLabels, ...
                                 flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric,'pdc',strID)

%%
%
% PDC inference p-values

flgPrinting  =   [1 1 1 2 3 0 0];
flgScale = 2;
[hxlabel hylabel] = xplot_pvalues(strID, c, ...
                                  flgPrinting,fs,w_max,chLabels,flgColor,flgScale);
xplot_title(alpha,metric,'p-value PDC',strTitle);

%% iDTF estimation
%
% metric = 'info';
% 

d = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
    d.Tragct = Tr_gct;
    d.pvaluesgct = pValue_gct;
%%
%
% |iDTF|^2 Matrix Layout Plotting
%

flgPrinting = [1 0 1 1 0 0 2]; % overriding default setting
flgScale = 1;
flgMax = 'TCI';
flgSignifColor = 3;
flgColor = [0];

w_max=fs/2;

strTitle1 = ['Model I with feedback: '];
strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p)  ' ]'];
strTitle =[strTitle1 strTitle2];

strID = 'Baccala & Sameshima (2001) Example 3';
[hxlabel,hylabel] = xplot(strID,d,flgPrinting,fs,w_max,chLabels, ...
                                 flgColor,flgScale,flgMax,flgSignifColor);

xplot_title(alpha,metric,'dtf',strTitle);

%%
%
% |iDTF|^2 Matrix Layout Plotting with Fine y-axis scaling
%

flgPrinting = [1 1 1 2 0 0 2]; % overriding default setting
flgScale = 2;
flgMax = 'TCI';
flgSignifColor = 3;
flgColor = [1];

w_max=fs/2;

strTitle1 = ['Model I with feedback: '];
strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p)  ' ]'];
strTitle =[strTitle1 strTitle2];

strID = 'Baccala & Sameshima (2001) Example 3';
[hxlabel,hylabel] = xplot(strID,d,flgPrinting,fs,w_max,chLabels, ...
                                 flgColor,flgScale,flgMax,flgSignifColor);

xplot_title(alpha,metric,'dtf',strTitle);


% End of example simulation.
