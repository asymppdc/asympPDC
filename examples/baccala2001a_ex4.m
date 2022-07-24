%% BACCALA & SAMESHIMA (2001A) EXAMPLE 4
% DESCRIPTION:
% Five-dimensional linear VAR[2] Model Example 4
%
%    x1==>x2  x2-->x3 x3-->x4 x4-->x5 x5-->x4
%
% Example taken from Baccala & Sameshima. Partial directed coherence: a new 
% concept in neural structure determination. 
% _Biol. Cybern._ *84*:463-474, 2001.
%
%                <http://dx.doi.org/10.1007/PL00007990>

%% Data sample generation

clear; clc; format compact; format short

nDiscard = 1000;    % number of points discarded at beginning of simulation
nPoints  = 2000;   % number of analyzed samples points

u = fbaccala2001a_ex4( nPoints, nDiscard );
chLabels = []; % or  = {'x_1';'x_2';'x_3';'x_4';'x_5'};
fs = 1;

%% Interaction diagram
%
% <<fig_baccala2001a_ex4_graph.png>>
%
% Figure 3a from Baccala & Sameshima. _Biol. Cybern._ *84*:463-474, 2001.

%% Equation
%
% <<fig_baccala2001a_ex4_eq.png>>
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
flgPrintResults = 1; % Flag to control printing gct_alg.m results on command window.
[Tr_gct, pValue_gct, Tr_igct, pValue_igct] = gct_alg(u,A,pf,gct_signif, ...
                                              igct_signif,flgPrintResults);
 
%% Original PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_dtf.m or issue 
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
flgColor = [0];
w_max=fs/2;
flgPrinting = [1 1 1 0 0 0 3]; % plot auto PDC on main diagonal

% h=figure;
% set(h,'NumberTitle','off','MenuBar','none', ...
%     'Name', 'Baccala & Sameshima (2001): Example 4')

[hxlabel hylabel] = xplot('Baccala & Sameshima (2001): Example 4',c,...
    flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha,metric,'pdc');

%% Original DTF estimation
%
% PDC analysis results are saved in *d* structure.
% See asymp_dtf.m or issue 
%
%   >> help asymp_dtf 
%
% for more detail.
metric = 'euc';
d = asymp_dtf(u,A,pf,nFreqs,metric,alpha);

%%
% $|DTF(\lambda)|^2$ Matrix Layout Plotting
flgColor = [0];
w_max=fs/2;

% h=figure;
% set(h,'NumberTitle','off','MenuBar','none', ...
%     'Name', 'Baccala & Sameshima (2001): Example 4')
vBarTitle = 'Baccala & Sameshima (2001): Example 4';
[hxlabel hylabel] = xplot(vBarTitle,d,...
    flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha,metric,'dtf',vBarTitle);

%%
% * Check & compare this results with Fig. 3b, page 469, Baccala & Sameshima (2001).
% * Note that in the article the amplitude,|PDC|, has been shown.
%    Here we preferred to plot $|PDC|^2$ and $|DTF|^2$.
%

%% Result from the original article, Baccala & Sameshima (2001) 
% In the original article the results is as following for comparison.
%
% 
% <<fig_baccala2001a_ex4ab.png>>
% 

%%
%