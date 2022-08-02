%% BACCALA & SAMESHIMA (2001A) EXAMPLE 4
%
% DESCRIPTION:
%
% Five-dimensional linear VAR[2] Model Example 4
%
%    $x1==>x2  x2-->x3 x3-->x4 x4<-->x5$
%
% Example from
%
%        Baccala & Sameshima (2001). Partial directed coherence: a new
%        concept in neural structure determination. _Biol. Cybern._
%        *84*:463-474.
%
%                <https://dx.doi.org/10.1007/PL00007990>
%% Other routines
%  See also  mvar, mvarresidue, asymp_pdc, asymp_dtf, gct_alg, igct_alg, 
%            xplot, xplot_pvalues
%  <baccala2001a_ex4.html |baccala2001a_ex4|> |
%%
clear; clc; format compact; format short

%% Interaction diagram
%
% <<fig_baccala2001a_ex4_graph.png>>
%
% Figure 3a from Baccala & Sameshima. _Biol. Cybern._ *84*:463-474, 2001.

%% Equation 5-dimension VAR[2] with feedback
%
% <<fig_baccala2001a_ex4_eq.png>>
%

%% Data sample generation

nDiscard = 1000;   % number of points to be discarded at beginning of simulation
nPoints  = 2000;   % number of analyzed samples points

u = fbaccala2001a_ex4( nPoints, nDiscard );
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
if flgDetrend
   for i=1:nChannels, u(i,:)=detrend(u(i,:)); end
   disp('Time series were detrended.');
end
if flgStandardize
   for i=1:nChannels, u(i,:)=u(i,:)/std(u(i,:)); end
   disp('Time series were scale-standardized.');
end

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
MVARadequacy_signif = 0.05; % VAR model estimation adequacy significance level
aValueMVAR = 1 - MVARadequacy_signif; % Confidence value for the testing

flgPrintResults = 1;
[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueMVAR,h,...
                                           flgPrintResults);

%% Granger causality test (GCT) and instantaneous GCT
%

gct_signif  = 0.01;  % Granger causality test significance level
igct_signif = 0.01;  % Instantaneous GCT significance level

metric = 'diag';     % euc  - original PDC or DTF;
                     % diag - generalized PDC (gPDC) or DC;
                     % info - information PDC (iPDC) or iDTF.
flgPrintResults = 1; % To print or not GCT/iGCT results on command window.

[Tr_gct, pValue_gct]   =  gct_alg(u,A,pf, gct_signif,flgPrintResults);
[Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgPrintResults);


%% Original PDC definition estimation
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

%% $|PDC(\lambda)|^2$  Matrix-Layout Plotting
%

flgColor = [0];
w_max=fs/2;
flgPrinting = [1 1 1 0 0 0 3]; % plot auto PDC on main diagonal

[h1,~, ~] = xplot('Baccala & Sameshima (2001a) - Example 4',c,...
    flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha,metric,'pdc');


%% Original DTF definition estimation
%
% DTF analysis results are saved in *d* structure.
% See asymp_dtf.m.

metric = 'euc';
d = asymp_dtf(u,A,pf,nFreqs,metric,alpha);

%% $|DTF(\lambda)|^2$  Matrix Layout Plotting
%

flgColor = [0];
w_max=fs/2;

vBarTitle = 'Baccala & Sameshima (2001) - Example 4';
[h2,~, ~] = xplot(vBarTitle,d,flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha,metric,'dtf',vBarTitle);

%%
% * Check & compare this results with Fig. 3b, page 469, Baccala & Sameshima (2001).
% * Note that in the article the amplitude, $|PDC|$, has been depicted.
%   While here we preferred to plot $|PDC|^2$ and $|DTF|^2$.
%

%% Result from the original article, Baccala & Sameshima (2001) 
% In the original article the results is as follow for comparison.
%
% 
% <<fig_baccala2001a_ex4ab.png>>
% 

%%
%