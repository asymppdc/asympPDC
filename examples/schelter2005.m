%% Schelter et al. (2005)  Five-dimension VAR[4]-process 
%
% Schelter, Winterhalder, Eichler, Peifer,Hellwig, Guschlbauer, L?cking,
% Dahlhaus & Timmer. Testing for directed influences among neural signals 
% using partial directed coherence. J. Neurosci Methods 152:210-9, 2005.
%
% <https://dx.doi.org/10.1016/j.jneumeth.2005.09.001>
% 
% Example Eq. (5) Five-dimension VAR[4]-process 
%


%% Data sample generation

clear; clc
flgPrintScreen = 'screen';

nDiscard = 5000;    % number of points discarded at beginning of simulation
nPoints  = 2000;   % number of analyzed samples points
N = nDiscard + nPoints; % number of simulated points

disp('======================================================================');
disp('       Schelter et al. J Neurosci Methods. 152:210-9, 2005.')
disp('               Linear 5-dimension VAR[4]-process')
disp('       x2==>x1  x3-->x2 x3==>x4 x3-->x5 x4==x2 x5-->x3  x5==x4');
disp('======================================================================');

randn('state', sum(100*clock))
% Variables initialization
ei=randn(5,N);
x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);
x4=zeros(1,N);
x5=zeros(1,N);
for t=1:4,
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1); x4(t)=randn(1);
   x5(t)=randn(1);
end;

%chLabels = []; % or 
chLabels = {'x_1';'x_2';'x_3';'x_4';'x_5'};

for t=5:N,
   x1(t) = 0.6*x1(t-1) + 0.65*x2(t-2) + ei(1,t);
   x2(t) = 0.5*x2(t-1) - 0.3*x2(t-2) - 0.3*x3(t-4) + 0.6*x4(t-1) + ei(2,t);
   x3(t) = 0.8*x3(t-1) - 0.7*x3(t-2) - 0.1*x5(t-3) + ei(3,t);
   x4(t) = 0.5*x4(t-1) + 0.9*x3(t-2) + 0.4*x5(t-2) + ei(4,t);
   x5(t) = 0.7*x5(t-1) - 0.5*x5(t-2) - 0.2*x3(t-1) + ei(5,t);
end;

y=[x1' x2' x3' x4' x5']; % data must be organized column-wise
u=y(nDiscard+1:N,:);

[nSegLength,nChannels]=size(u);

fs = 1;

%%
% Data pre-processing: detrending and standardization options

flgDetrend = 1;     % 1: Detrending the data set
flgStandardize = 0; % 0: No standardization
% Checking data dimension
[nChannels,nSegLength] = size(u);
if nChannels > nSegLength, 
   u = u.'; 
   [nChannels,nSegLength] = size(u);
end;
% Detrending 
if flgDetrend,
   for i = 1:nChannels, u(i,:) = detrend(u(i,:)); end;
   disp('Time series were detrended.');
end;
% Standardization 
if flgStandardize,
   for i = 1:nChannels, u(i,:) = u(i,:)/std(u(i,:)); end;
   disp('Time series were scale-standardized.');
end;

%%
% MVAR model estimation

maxIP = 30;         % maximum model order to consider.
alg = 1;            % 1 = Nutall-Strand MVAR estimation algorithm
criterion = 1;      % 1 = AIC, Akaike Information Criteria

disp('Running MVAR estimation and GCT analysis routines.')

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);

disp(['Number of channels = ' int2str(nChannels) ' with ' ...
    int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%%
%Testing the adequacy of MAR model fitting through Portmanteau test
h = 20; % testing lag value
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
                                                       
%% PDC, threshold and confidence interval calculations.
%
metric = 'diag'; % Generalized PDC estimation
nFreqs = 128;   % Number of frequency points to consider
alpha  = 0.0001;  % Significance level for PDC/DTF null hypothesis test

%%
% PDC analysis results are saved in *c* data structure.
% See asymp_dtf.m or issue 
%
%   >> help asymp_pdc 
%
% command for more detail.

c = asymp_pdc(u,A,pf,nFreqs,metric,alpha);

%% 
% Plotting options set up, mostly cosmetics, used in xplot.m routine:
%
switch lower(flgPrintScreen)
   case 'print'
      flgMax = 'TCI';
      flgSignifColor = 1; % black + gray
      flgScale = 3;       % [0 max(flgMax)]
   otherwise % e.g., 'screen'
      flgMax = 'TCI';
      flgSignifColor = 3; % red + green
      flgScale = 2;       % [0 1]/[0 .1]/[0 .01]
end;
flgScale
%%
% flgColor parameter for PDC matrix-layout plot.
flgColor = [0];   % Plotting option for automatic scaling for small PDC
                  % values.
                  % if flgColor = 0, y-axis scale = [0 1]
                  % elseif flgColor = 1, xplot routine rescale 
                  % the y-axis automatically according to following rules:
                  %   If .001<= max(|PDC(f)|^2) < .01 background-color = light-blue,
                  %                          so that y-axis scale = [0 .1]
                  %   elseif max(|PDC(f)|^2) < .001 background-color = light-purple
                  %                          and y-axis = [0 .01].

%             [1 2 3 4 5 6 7]
flgPrinting = [1 1 1 2 2 0 5];
%    blue-line | | | | | | 7--Spectra (0: w/o; 1: Linear; 2: Log; 3: PDC2; 
%              | | | | | |      4: Linear normalized; 5: Log spectra + PDC2)
%         gray | | | | | 6--Coh2 (0: w/o Coh2; 1: w Coh2)
%  dashed-blue | | | | 5--Plot lower confidence limit (***legacy)
%  dashed-blue | | | 4--Plot confidence interval
%         red  | | 3--Significant PDC2|DTF2 in red lines (***legacy)
% dashed-black | 2--Patnaik threshold level in black dashed-lines
%        green 1-- PDC2/DTF2 in green lines or black w/o statistics,
%                  See flgSignifColor for line color selection.

w_max=fs/2;

for kflgColor = flgColor,
   %    h=figure;
   %    set(h,'NumberTitle','off','MenuBar','none', ...
   %       'Name', 'Schelter et al. J. Neurosci Methods (2005)')
   %    [ax,hT]=suplabel(['Linear pentavariate, VAR[4]-process: ' ...
   %       int2str(nPoints) ' data points.'],'t');
   %    set(hT,'FontSize',12); % Title font size
   %   [hxlabel hylabel] = xplot(h,c,...
   %                               flgPrinting,fs,w_max,chLabels,kflgColor);
   
   strID = 'Schelter et al. J. Neurosci Methods (2005)';
   [hxlabel,hylabel] = xplot(strID,c,flgPrinting,fs,w_max,chLabels, ...
                                      kflgColor,flgScale,flgMax,flgSignifColor);
   xplot_title(alpha,metric,'pdc',strID);
end;

%% Remarks:
% 
% * Check the plot with Figs. 1 and 3 (see pages 212 and 215 Schelter
%   et al., 2005). As you may notice, the power spectrum of x4, as well 
%   as |PDC|2_24, differ significantly.
%   Our guess is that Schelter et al. may have used slightly different
%   parameters from what they stated for Eq. 5.
%
% * Note that, for linear model with balanced innovation, the maximum
%   of PDC estimates is roughly proportional to the autoregressive 
%   model coefficients. As we are plotting |PDC|^2, amplitude is 
%   roughly proportional to the square of VAR coefficients.
