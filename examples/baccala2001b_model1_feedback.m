%% Baccala & Sameshima (2001b) 7-dimensional VAR[2] model with loop and feedback
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
% # Koichi Sameshim, Daniel Y. Takahashi, Luiz A. Baccala. On the
% statistical performance of Granger-causal connectivity estimators. Brain
% Informatics (2015) 2:119?133.
%
% <http://dx.doi.org/10.1007/s40708-015-0015-1>
% 
% Example Model 1 - 7-dimensional VAR[2] model with loop and feedback

%% Data sample generation

clear; clc; format compact; format short

nDiscard = 5000;    % number of points discarded at beginning of simulation
nPoints  = 2000;   % number of analyzed samples points
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

metric = 'diag'; % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
flgPrintResults = 1;

% [Tr_gct, pValue_gct, Tr_igct, pValue_igct] = gct_alg(u,A,pf,gct_signif, ...
%                                               igct_signif,flgPrintResults);
[Tr_gct, pValue_gct]   = gct_alg (u,A,pf,gct_signif, flgPrintResults);
[Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgPrintResults);

%% Original PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_pdc.m or issue 
%
%   >> help asymp_pdc 
%
% command for more detail.
nFreqs = 128;
metric = 'euc';
alpha = 0.01;
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics

%%
% $|PDC(\lambda)|^2 Matrix Layout Plotting

flgPrinting = [1 1 1 2 2 1 2]; % overriding default setting
flgColor = 0;
w_max=fs/2;

strTitle1 = ['7-dimensional linear VAR[3] Model I with feedback: '];
strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p) '; ' ...
   datestr(now) ']'];
strTitle =[strTitle1 strTitle2];
% h=figure;
% set(h,'NumberTitle','off','MenuBar','none', ...
%    'Name', 'Baccala & Sameshima (2001) Example 3')

strTask = 'Baccala & Sameshima (2001) Example 3';
[hxlabel hylabel] = xplot(strTask,c,...
                          flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha,metric,strTask);
% [ax,hT]=suplabel( strTitle, 't' );
% set(hT,'FontSize',10)

%%
% The spectral coherence are plotted in gray-line. You may notice that
% although the isolated structures 6 and 7 have power peak at the same
% frequency of remaining structures, they most likely will present low or no
% coherence with other structures. The red-lines indicate statistically significant PDC. You may
% also see occasionally spurious false positive connectivity inference
% which may occur in approximately equal to the level of significance
% chosen for the null hypothesis for non causality test, i.e. $\alpha$.

%% 
% Theoretical PDC results from the original article, Baccala & Sameshima (2001b) 
% Figure 6b, from article.
%
% <<fig_baccala2001b_pdc_result.png>>
% 

%%
% In this original article's figure the significant PDC is shown in shaded
% area, and the spectral coherence in gray-line.

%% Generalized PDC estimation
%
% PDC analysis results are saved in *d* structure.
% See asymp_pdc.m or issue
%
%   >> help asymp_pdc 
%
% command to see for more detail.
nFreqs = 128;
metric = 'diag';
alpha = 0.01;
d = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics

%%
% PDCn Matrix Layout Plotting

flgPrinting = [1 1 1 2 2 0 2];
flgColor = [0];
w_max=fs/2;
flgScale = 1; % y-axis = [0 1]
flgMax = 'TCI';
flgSignifColor = 3; % red = significant, gree = nonsignificant
for kflgColor = flgColor,
%    h=figure;
%    set(h,'NumberTitle','off','MenuBar','none', ...
%       'Name', 'Baccala & Sameshima (2001b) Linear Model I')
   [hxlabel,hylabel] = xplot(strTask,d,flgPrinting,fs,w_max,chLabels, ...
                               kflgColor,flgScale,flgMax,flgSignifColor);
   xplot_title(alpha,metric,strTask);
%    [ax,hT]=suplabel(['Linear model ' ...
%       int2str(nPoints) ' data points.'],'t');
%    set(hT,'FontSize',10); % Title font size
end;


%% Generalized DTF = DC estimation
%
% DC analysis results are saved in *e* structure.
% See asymp_dtf.m or issue
%
%   >> help asymp_dtf 
%
% command for more detail.
nFreqs = 128;
metric = 'diag';
alpha = 0.01;
e = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics

%%
% DTF Matrix Layout Plotting

flgPrinting = [1 1 1 2 2 0 2];
flgColor = [1];
w_max=fs/2;
flgScale = 2; % y-axis = [0 1]
flgMax = 'TCI';
flgSignifColor = 3; % red = significant, gree = nonsignificant
for kflgColor = flgColor,
%    h=figure;
%    set(h,'NumberTitle','off','MenuBar','none', ...
%       'Name', 'Baccala & Sameshima (2001b) Linear Model I with feedback')
   [hxlabel,hylabel] = xplot(strTask,e,flgPrinting,fs,w_max,chLabels, ...
                               kflgColor,flgScale,flgMax,flgSignifColor);
   xplot_title(alpha,metric,'dtf',strTask);
%    [ax,hT]=suplabel(['Linear model ' ...
%       int2str(nPoints) ' data points.'],'t');
%    set(hT,'FontSize',10); % Title font size
end;
%%
% In this example with feedback, theoretically, all structures {1,2,3,4,5} 
% should be able to reach each other, henceforth all DC between these set 
% of structures should be significant. But false negative DC connectivity
% may occur more often for short simulation data segment.

%%
% * In the original article the amplitude of PDC/DTF has been plotted.
% Here we chose squared PDC, i.e. the xplot.m routine can only handle
% squared magnitude measures, $|PDC|^2$, $|DTF|^2$, $|Coh|^2$, as the 
% asymptotic statistics is formulated in spectral coherence domain. 

