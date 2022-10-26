%% BACCALA & SAMESHIMA (2001a) - Example 3, p.468
%
% DESCRIPTION:
%
% Linear five-dimension  VAR(3) model with feedback
%
% Baccala & Sameshima. Partial directed coherence: a new concept in neural
% structure determination. _Biol. Cybern._ *84*:463-474, 2001.
%
% <http://dx.doi.org/10.1007/PL00007990>
%                
% Example 3 (pag.468)
%                       VAR(3) with feedback between x4 and x5
%
%     $x1-->x2  x1-->x3 x1-->x4 x4<-->x5$
%
% Example *baccala2001_ex3.m*:
%
%% See also: mvar, mvarresidue, asymp_pdc, asymp_dtf, gct_alg, 
%              igct_alg, xplot, xplot_pvalues             

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.

%%

clear; clc; format compact; format short

%% Interaction diagram
%
% <<fig_baccala2001a_ex3_graph.png>>
%
% Figure 2a from Baccala & Sameshima. _Biol. Cybern._ *84*:463-474, 2001.

%% Equation system
%
% <<fig_baccala2001a_ex3_eq.png>>
%

%% Data sample generation
%
nBurnIn = 5000;   % number of points discarded at beginning of simulation
nPoints  = 1000;   % number of analyzed samples points
u = fbaccala2001a_ex3(nPoints, nBurnIn); % Model function
chLabels = []; %{'x_1';'x_2';'x_3';'x_4';'x_5'};
fs = 1; % Normalized frequency.


%% Data pre-processing: detrending and/or normalization options
%

[nChannels,nSegLength] =size(u);
if nChannels > nSegLength
   u = u.'; 
   [nChannels,nSegLength]=size(u);
end

for i=1:nChannels, u(i,:)=detrend(u(i,:)); end
disp('Time series were detrended.');

%% MVAR model estimation

maxIP = 30;         % maximum model order to be considered.
alg = 1;            % 1: Nutall-Strand MVAR estimation algorithm
criterion = 1;      % 1: AIC, Akaike Information Criteria

disp('Running MVAR estimation routine.')

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);

disp(['Number of channels = ' int2str(nChannels) ' with ' ...
    int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%%
% Testing for adequacy of MAR model fitting through Portmanteau test

h = 20; % GCT testing lag
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

flgPrintResults = 1;

[Tr_gct, pValue_gct] = gct_alg(u,A,pf,gct_signif,flgPrintResults);
[Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgPrintResults);

                                           
%% Original PDC estimation
%
% PDC analysis results are saved in *c* struct variable.
% See asymp_pdc.m or issue 
%   >> help asymp_pdc 
% to see more detail.

nFreqs = 128;
metric = 'euc'; % Euclidean metric = original PDC definition
alpha = 0.01;
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha);

c.Tragct = Tr_gct;         % Assigning GCT results to c struct variable.
c.pvaluesgct = pValue_gct;
    
%% $|PDC(\lambda)|^2$ Matrix-Layout Plot

flgPrinting = [1 1 1 2 3 0 3]; % overriding default setting
flgColor = 0;
w_max=fs/2;

strTitle = ['Linear 5-dim VAR(3) model w feedback'];  %'; '  datestr(now) ']'];
strBarTitle = 'Baccala & Sameshima (2001a) Example 3';

[h1,~, ~] = xplot(strBarTitle,c,...
                          flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha,metric,'pdc', strTitle);

%%
% Original PDC2 plots with auto-PDC2 along the main-diagonal, reproducing
% the figure of original article. Note however that here we are depicting the
% squared-PDC, while in the article |PDC| was used.


%% information PDC estimation
%
% iPDC analysis results are saved in *d* struct variable.
% See asymp_pdc.m 

nFreqs = 128;
metric = 'info'; % Choosing information PDC

d = asymp_pdc(u,A,pf,nFreqs,metric,alpha);

d.Tragct = Tr_gct;
d.pvaluesgct = pValue_gct;


%% $|_{i}PDC(\lambda)|^2$ Matrix-Layout Plotting

%                 1 2 3 4 5 6 7
   flgPrinting = [1 1 1 2 3 0 2]; % GCT and power spectra selection.
%                 | | | | |   |     
%           blue  | | | | |   7-- 2: Log spoectra;
%    dark-purple  | | | | 5-- Print GCT p-values and dot-mark significant
%  or dark-green  | | | |     connectivity channel-pair 3: p-values + dot-mark 
%                 | | | |     significant GCT
%    dashed-blue  | | | 4-- Confidence interval 2: shaded-plot
%            red  | | 3-- Significant PDC2 3: in red lines
%   dashed-black  | 2-- Patnaik threshold level in black dashed-lines
%          green  1-- PDC2 in green lines or black w/o statistics,
%                           see flgSignifColor bellow for line color selection.

flgColor = 0; w_max=fs/2;

[h2,~, ~] = xplot(strBarTitle,d,...
                          flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha,metric,'pdc', strTitle);

%%
% Information PDC2 plots with log-spectra along the main-diagonal.

%% Original DTF estimation
%
% DTF analysis results are saved in *e* struct variable.  
% See asymp_dtf.m for more detail.

metric = 'euc'; % Original DTF
e = asymp_dtf(u,A,pf,nFreqs,metric,alpha);


%% $|{DTF}(\lambda)|^2$ Matrix-Layout Plot
%

flgColor = 0; flgScale = 1; flgMax = 'dtf'; flgSignifColor = 3;
flgPrinting = [1 1 1 2 0 0 3]; % overriding default setting

[h3,~,~] = xplot(strBarTitle,e,flgPrinting,fs,w_max,chLabels, ...
                                 flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric,'dtf',strTitle);

%%
% Original DTF2 plots with auto-DTF2 along the main-diagonal,
% closely reproducing the figure of original article in which we depicted |DTF|. 


%% Result from the original article, Baccala & Sameshima (2001) 
% Figure 2, page 469 from article.
%
% <<fig_baccala2001a_ex3abc.png>>
% 

%% Some remarks:
% 
% # Check & compare Fig. 2b, page 468, Baccala & Sameshima (2001).
% # Note that in the original article the amplitude PDC has been plotted,
%   while here we preferred to plot squared-PDC.

%%
% This completes the Example 3 (Baccala & Sameshima, 2001)
