%% BACCALA & SAMESHIMA (2001a) - Example 5
%
% DESCRIPTION:
%
% Five-dimensional linear VAR[2] Model Closed-loop Example 5, pages 468 and 469,
%
%    x1==>x2  x2-->x3 x3-->x4 x4<-->x5 x5-->x1
%
% borrowed from:
%       Baccala & Sameshima (2001). Partial directed coherence: a new concept in neural
%       structure determination. _Biol. Cybern._ *84*:463--474.
%
% <http://dx.doi.org/10.1007/PL00007990>
% 
% Example Five-dimensional VAR[2] with loop and feedback
%
%% Other routines
%  See also  mvar, mvarresidue, asymp_pdc, asymp_dtf, gct_alg, igct_alg, 
%            xplot, xplot_pvalues
%  <baccala2001a_ex5.html |baccala2001a_ex5|> |
%%

clear; clc; format compact; format short


%% Interaction diagram
%
% <<fig_baccala2001b_example5_graph.png>>
%
% Figure 4a from Baccala & Sameshima. _Biol. Cybern._ *84*:463-474, 2001.

%% Equation 5-dimension VAR[2] with closed-loop and feedback
%
% <<fig_baccala2001b_example5_eq.png>>
%

%% Data sample generation
% 
nDiscard = 10000;   % number of points discarded at beginning of simulation
nPoints  = 5000;    % number of analyzed samples points
u = fbaccala2001a_ex5( nPoints, nDiscard );
chLabels = []; %{'x_1';'x_2';'x_3';'x_4';'x_5'};
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

% Detrending the data set and no detrending
for i=1:nChannels, u(i,:)=detrend(u(i,:)); end
disp('Time series were detrended,');
disp('                 and not scale-standardized.');

%% MVAR model estimation
%

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

%% Granger causality test (GCT) and instantaneous GCT
%

gct_signif  = 0.01;  % Granger causality test significance level
igct_signif = 0.01;  % Instantaneous GCT significance level
flgPrintResults = 1; % Flag to control printing gct_alg.m results on command window.

[Tr_gct, pValue_gct]   =  gct_alg(u,A,pf, gct_signif, flgPrintResults);
[Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgPrintResults);

%% Original PDC estimation
% PDC analysis results are saved in *c* struct variable.
% See asymp_pdc.m for details.

nFreqs = 256;
metric = 'euc';  % Euclidean = original PDC (PDC) or DTF.
alpha = 0.01;

c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
c.pvaluesgct = pValue_gct; % Necessary as flgPrinting(5) = 3, i.e. printing GCT
c.Tragct = Tr_gct;

%% iPDC2 Matrix-Layout Plotting

flgPrinting = [1 1 1 0 3 0 3]; % With GCT and log-spectra on main diagonal
flgColor = 0;
w_max=fs/2;

strBarTitle = 'Baccala & Sameshima (2001A) - Example 5';
strTitle = ['Linear model closed-loop: ' int2str(nPoints) ' data points.']

[h1,hxlabel hylabel] = xplot(strBarTitle,c,flgPrinting,fs,w_max,chLabels,flgColor);

xplot_title(alpha,metric,'pdc',strTitle);
pause(5)

%% Original DTF estimation
% DTF analysis results will be saved in *d* struct variable.
% See asymp_dtf.m for further details.

metric = 'euc';
d = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Estimate DTF and asymptotic statistics


%% DTF2 Matrix Layout Plotting with fixed y-axis scale
flgPrinting = [1 1 1 0 0 0 3]; % Plot log-spectra on main-diagonal
flgColor = 1;
w_max=fs/2;
flgMax = 'TCI';
flgScale = 1;
flgSignifColor = 3;

[h2,hxlabel,hylabel] = xplot(strBarTitle,d,flgPrinting,fs,w_max,chLabels, ...
                                       flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric,'dtf',strTitle);
pause(5)

%%
% Note that the magnitude is not necessarily adequate criteria for the
% presence of connectivity. 

%% DTF Matrix Layout Plotting with y-axis scaling
%

flgPrinting = [1 1 1 2 2 0 0]; % Without power spectra on main-diagonal
flgColor = 1;
w_max=fs/2;
flgMax = 'TCI';
flgScale = 3;
flgSignifColor = 3;  

[h3,hxlabel,hylabel] = xplot(strBarTitle,d,flgPrinting,fs,w_max,chLabels, ...
                               flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric,'dtf',strTitle);
pause(5)
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
