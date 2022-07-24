%% Baccala & Sameshima (2001b) 7-dimension VAR[2] model with loop and feedback
%
% Description:
%
% Baccala & Sameshima. Overcoming the limitations of correlation analysis 
% for many simultaneously processed neural structures, Progress in Brain 
% Research, 130:33--47, 2001.
%
% <http://dx.doi.org/10.1016/S0079-6123(01)30004-3>
%
% See also:
%
% 
% # Koichi Sameshima ? Daniel Y. Takahashi ? Luiz A. Baccala ?. On the
% statistical performance of Granger-causal connectivity estimators. Brain
% Informatics (2015) 2:119?133.
%
% <http://dx.doi.org/10.1007/s40708-015-0015-1>
% 
% Example Model 1 - 7-dimension VAR[2] model with loop and feedback

%% Data sample generation

clear; clc; format compact; format short

nDiscard = 50000;    % number of points discarded at beginning of simulation
nPoints  = 10000;   % number of analyzed samples points
N=nDiscard+nPoints; % number of simulated points

u = fbaccala2001b_model1_feedback( nPoints, nDiscard );

[nSegLength,nChannels]=size(u);

%chLabels = []; % or 
chLabels = {'x_1';'x_2';'x_3';'x_4';'x_5';'x_6';'x_7'};

fs = 1;
%% Interaction diagram
%
% <<fig_baccala2001b_graph.png>>
%
% Figure 2a from Baccala & Sameshima. _Biol. Cybern._ *84*:463-474, 2001.

%% Equation Model I with feedback
%
% <<fig_baccala2001b_eq.png>>
%

%%
% Data pre-processing: detrending and normalization options

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

alpha = 0.001;

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

%%
% Granger causality test (GCT) and instantaneous GCT

gct_signif  = 0.01;  % Granger causality test significance level
igct_signif = 0.01;  % Instantaneous GCT significance level

metric = 'info'; % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
flgPrintResults = 1;

[Tr_gct, pValue_gct, Tr_igct, pValue_igct] = gct_alg(u,A,pf,gct_signif, ...
                                              igct_signif,flgPrintResults);
                                                       
%% Original PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_pdc.m or issue 
%
%   >> help asymp_pdc 
%
% command for more detail.
nFreqs = 128;
% 
% c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
% 
% %%
% % $|_iPDC(\lambda)|^2 Matrix Layout Plotting
% 
% flgPrinting = [1 1 1 2 2 0 0]; % overriding default setting
% flgScale = 1;
% flgMax = 'pdc';
% flgSignifColor = 3;
% flgColor = [0];
% 
% w_max=fs/2;
% 
% strTitle1 = ['Model I with feedback: '];
% strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p)  ' ]'];
% strTitle =[strTitle1 strTitle2];
% strID = 'Baccala & Sameshima (2001) Example 3';
% [hxlabel,hylabel] = xplot(strID,c,flgPrinting,fs,w_max,chLabels, ...
%                                  flgColor,flgScale,flgMax,flgSignifColor);
% xplot_title(alpha,metric,'pdc',strTitle);

% % [ax,hT]=suplabel( strTitle, 't' );
% % set(hT,'FontSize',16)



metric = 'info';

d = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics

%%
% $|_iPDC(\lambda)|^2 Matrix Layout Plotting

flgPrinting = [1 1 1 2 2 0 0]; % overriding default setting
flgScale = 1;
flgMax = 'd';
flgSignifColor = 3;
flgColor = [0];

w_max=fs/2;

strTitle1 = ['Model I with feedback: '];
strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(d.p)  ' ]'];
strTitle =[strTitle1 strTitle2];

% h=figure;
% set(h,'NumberTitle','off','MenuBar','none', ...
%    'Name', 'Baccala & Sameshima (2001) Example 3')

strID = 'Baccala & Sameshima (2001) Example 3';
[hxlabel,hylabel] = xplot(strID,d,flgPrinting,fs,w_max,chLabels, ...
                                 flgColor,flgScale,flgMax,flgSignifColor);

xplot_title(alpha,metric,'dtf',strTitle);
% [ax,hT]=suplabel( strTitle, 't' );
% set(hT,'FontSize',16)