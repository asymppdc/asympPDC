%% BACCALA & SAMESHIMA (2001A) EXAMPLE 5
% DESCRIPTION:
%
% Five-dimensional linear VAR[2] Model Example 4
%
%    x1==>x2  x2-->x3 x3-->x4 x4<-->x5 x5-->x1
%
% Example borrowed from:
%       Baccala & Sameshima (2001). Partial directed coherence: a new concept in neural
%       structure determination. _Biol. Cybern._ *84*:463--474.
%
% <http://dx.doi.org/10.1007/PL00007990>
% 
% Example Five-dimensional VAR[2] with loop and feedback

close all; clear; clc

%% Data sample generation
% 
nDiscard = 10000;    % number of points discarded at beginning of simulation
nPoints  = 2000;    % number of analyzed samples points
u = fbaccala2001a_ex5( nPoints, nDiscard );
chLabels = []; % or  = {'x_1';'x_2';'x_3';'x_4';'x_5'};
fs = 1; 

%%
% Data pre-processing: detrending and normalization options

flgDetrend = 1;     % Detrending the data set
flgStandardize = 0; % No standardization
[nChannels,nSegLength] =size(u);
if nChannels > nSegLength
   u = u.'; 
   [nChannels,nSegLength]=size(u);
end

% Detrending the data set and no detrending
for i=1:nChannels, u(i,:)=detrend(u(i,:)); end
disp('Time series were detrended,');
disp('                 and not scale-standardized.');

%%
% MVAR model estimation

maxIP = 30;         % maximum model order to consider.
alg = 1;            % 1: Nutall-Strand MVAR estimation algorithm;
%                   % 2: minimum least squares methods;
%                   % 3: Vieira Morf algorithm;
%                   % 4: QR ARfit algorith.

criterion = 1;      % Criterion for order choice:
%                   % 1: AIC, Akaike Information Criteria; 
%                   % 2: Hanna-Quinn;
%                   % 3: Schwartz;
%                   % 4: FPE;
%                   % 5: fixed order given by maxIP value.

disp('Running MVAR estimation routine...')
[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);
pause(3);
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
flgPrintResults = 1; % Flag to control printing gct_alg.m results on command window.
% [Tr_gct, pValue_gct, Tr_igct, pValue_igct] = gct_alg(u,A,pf,gct_signif, ...
%                                               igct_signif,flgPrintResults);
[Tr_gct, pValue_gct] = gct_alg(u,A,pf,gct_signif,flgPrintResults);
[Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgPrintResults);

%% Information PDC estimation
% PDC analysis results are saved in *c* struct variable.
% See asymp_pdc.m for details.

nFreqs = 256;
metric = 'info';  % info = information PDC (iPDC) or iDTF.
alpha = 0.01;
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
c.pvaluesgct = pValue_gct; % Necessary as flgPrinting(5) = 3, i.e. printing GCT
c.Tragct = Tr_gct;

%%
% PDCn Matrix Layout Plotting

flgPrinting = [1 1 1 0 3 0 2]; % With GCT and log-spectra on main diagonal
flgColor = 0;
w_max=fs/2;
vTitle = 'Baccala & Sameshima (2001A) - Example 5';
[hxlabel hylabel] = xplot(vTitle,c,flgPrinting,fs,w_max,chLabels,flgColor);

xplot_title(alpha,metric,'pdc',['Linear model Example 5: ' ...
                                             int2str(nPoints) ' data points.']);

%% Original/Information DTF estimation
%
% DTF analysis results will be saved in *d* struct variable.
% See asymp_dtf.m for further details.

metric = 'info';
d = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Estimate DTF and asymptotic statistics


%%
% DTF Matrix Layout Plotting with fixed y-axis scale
flgPrinting = [1 1 1 0 0 0 2]; % Plot log-spectra on main-diagonal
flgColor = 1;
w_max=fs/2;
flgMax = 'TCI';
flgScale = 1;
flgSignifColor = 3;

[hxlabel,hylabel] = xplot(vTitle,d,flgPrinting,fs,w_max,chLabels, ...
                                       flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric,'dtf',['Linear model Example 5: ' ...
                                             int2str(nPoints) ' data points.']);

%%
% Note that the magnitude is not necessarily adequate criteria for the
% presence of connectivity. 

%%
% DTF Matrix Layout Plotting with y-axis scaling

flgPrinting = [1 1 1 2 2 0 0]; % Without power spectra on main-diagonal
flgColor = 1;
w_max=fs/2;
flgMax = 'TCI';
flgScale = 3;
flgSignifColor = 3;  

[hxlabel,hylabel] = xplot(vTitle,d,flgPrinting,fs,w_max,chLabels, ...
                               flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric,'dtf',['Linear model Example 5: ' ...
                                            int2str(nPoints) ' data points.']);

%%
% Note that theoretically any node or strutucture is reachable from all other
% structures. If sufficiently large sample size is used, all DTF will be
% significant.
%
% _Suggestion_: try  playing with different sample sizes, i.e. changing
% *nPoints* parameter.

%% Concluding remarks 
% * Check & compare with Fig.4b, page 469 in Baccala & Sameshima (2001).
% * In the original article the amplitude PDC has been plotted.
%   Here we preferred to graph squared-PDC and DTF.
