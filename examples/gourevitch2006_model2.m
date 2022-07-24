%% GOUREVITCH ET AL. (2006) LINEAR 2-VAR WITH FEEDBACK AND COMMON SOURCE
% DESCRIPTION:
%
% Example Model 2: Linear bivariate model with bidirectional influence 
%                          with common source
%
%     *x1<==>x2* (feedback)  *x1==x2* (instantaneous causality)
%
%     *x1 <== _S_ ==>x2*, which gives x1==x2 (iGC) from hidden source
%
% Example *Model 2* borrowed from:
% *Gourevitch, Bouquin-Jeannes & Faucon*. Linear and nonlinear casuality between 
%     signals: methods, examples and neurophysiological applications. 
%     _Biol Cybern_ *95*:349--369, 2006.
%
% <http://dx.doi.org/10.1007/s00422-006-0098-0>
%
%% Other routines
%  See also  asymp_pdc, asymp_dtf, gct_alg2, xplot
%            xplot_pvalues
% < gourevitch2006_model2.html |gourevitch2006_model2|> |

%%

clear; clc

%% Data sample generation
% 
nDiscard = 1000;    % number of points discarded at beginning of simulation
nPoints  = 800;     % number of analyzed samples points

u = fgourevitch2006_model2( nPoints, nDiscard );

chLabels = []; % or  = {'x_1';'x_2';'x_3';'x_4';'x_5'};
fs = 1; % normalized frequency

%% Interaction diagram
%
% <<fig_gourevitch2006_model2_graph.png>>
%
% Figure 4 from Gourevitch. _Biol. Cybern._ *95*:349-369, 2006.

%% Equation
%
% <<fig_gourevitch2006_model2_eq.png>>
%

%%
% Equations (38) and (39), Gourevitch et al. (2006)


%%
% Data pre-processing: detrending and normalization options

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
alg = 1;            % 1: Nutall-Strand MVAR estimation algorithm;
%                   % 2: minimum least squares methods;
%                   % 3: Vieira Morf algorithm;
%                   % 4: QR ARfit algorith.

criterion = 1;      % Criterion for order choice:
%                   % 1: AIC, Akaike Information Criteria; 
%                   % 2: Hanna-Quinn;
%                   % 3: Schwarz;
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
alpha = 0.01;

% Granger causality test (GCT) and instantaneous GCT

gct_signif  = alpha;  % Granger causality test significance level
igct_signif = alpha;  % Instantaneous GCT significance level
flgPrintResults = 1; % Flag to control printing gct_alg.m results on command window.
[Tr_gct, pValue_gct] = gct_alg2(u,A,pf,gct_signif,igct_signif,flgPrintResults);
 
%% Original PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_pdc.m or issue 
%   >> help asymp_pdc 
% command for more detail.

nFreqs = 128;
metric = 'euc';  % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
c.Tragct = Tr_gct;  c.pvaluesgct = pValue_gct;
%%
% PDC Matrix Layout Plotting

flgPrinting = [1 1 1 2 2 0 1]; % overriding default setting
flgColor = 1;
w_max=fs/2;
chLabels={'X_1';'X_2'}; %Optional channel labels;


for kflgColor = flgColor,
   %    h=figure;
   %    set(h,'NumberTitle','off','MenuBar','none', ...
   %          'Name', 'Gourevitch et al.(Biol Cybern, 2006)')
   
   strID = 'Gourevitch et al.(Biol Cybern, 2006)';
   [hxlabel hylabel] = xplot(strID,c,flgPrinting,fs,w_max,chLabels,kflgColor);
   strTitle = ['Model 2: Linear bivariate model with common source: ' ...
               int2str(nPoints) ' data points.'];
   xplot_title(alpha,metric,'pdc',strTitle); %Main title with PDC type and alpha value
end;

%% Generalized PDC estimation
%
% gPDC analysis results are saved in *d* structure.
% See asymp_dtf.m or issue 
%   >> help asymp_pdc 
% command for more detail.

metric = 'diag';  % euc  = original PDC or DTF;
                  % diag = generalized PDC (gPDC) or DC;
                  % info = information PDC (iPDC) or iDTF.
d = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
d.Tragct = Tr_gct;  d.pvaluesgct = pValue_gct;
%%
% PDC Matrix Layout Plotting

flgPrinting = [1 1 1 2 2 0 3]; % overriding default setting
flgColor = 1;
w_max=fs/2;
chLabels={'X_1';'X_2'}; %Optional channel labels;

   
for kflgColor = flgColor,
   %    h=figure;
   %    set(h,'NumberTitle','off','MenuBar','none', ...
   %          'Name', 'Gourevitch et al.(Biol Cybern, 2006)')
   
   strID = 'Gourevitch et al.(Biol Cybern, 2006)';
   [hxlabel hylabel] = xplot(strID,d,flgPrinting,fs,w_max,chLabels,kflgColor);
   %
   strTitle = ['Model 2: Linear bivariate model with common source: ' ...
                                              int2str(nPoints) ' data points.'];
   xplot_title(alpha,metric,'pdc',strTitle); %Main title with PDC type and alpha value
end;

%% Result from the original article, Gourevitch (2006) 
%
% 
% <<fig_gourevitch2006_model2_result.png>>
% 

%%
% *Figure* - from Gourevitch et al. (2006) presented results of PDC and DC


%% gPDC p-values matrix layout plots
%

flgPrinting  =   [1 1 1 2 3 0 0];
flgScale = 2;
[hxlabel hylabel] = xplot_pvalues([],d,flgPrinting,fs,w_max,chLabels, ...
                                                             flgColor,flgScale);

xplot_title(alpha,metric,['p-value gPDC'],'Gourevitch et al.(Biol Cybern, 2006)');

%%
% gPDC's p-values in the frequency domain, which, in this case, are < 1.0e-15.

%% gDTF = DC estimation
%
% DC analysis results will be saved in *e* structure.
% See asymp_dtf.m or issue 
%   >> help asymp_dtf 
% command for more detail.
metric = 'diag';
e = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Estimate DTF and asymptotic statistics

%%
% DTF Matrix Layout Plotting option with fixed y-axis scale on [0 1] range
%
flgPrinting = [1 1 1 2 2 0 3]; % Plot auto-DTF on main-diagonal
flgColor = 1; w_max=fs/2; flgMax = 'TCI'; flgScale = 1; flgSignifColor = 1;  

w_max=fs/2;
chLabels={'X_1';'X_2'}; %Optional channel labels;
for kflgColor = flgColor,
%    h=figure;
%    set(h,'NumberTitle','off','MenuBar','none', ...
%          'Name', 'Gourevitch et al.(Biol Cybern, 2006)')
   strID = 'Gourevitch et al.(Biol Cybern, 2006)';
   [hxlabel hylabel] = xplot(strID,e,flgPrinting,fs,w_max,chLabels,kflgColor);
   strTitle = ['Model 2: Linear bivariate model with common source: ' ...
               int2str(nPoints) ' data points.'];
   xplot_title(alpha,metric,'dtf',strTitle);
end;



%%
% Note significant Instantaneous Granger causality in between the pair 
%  of processes. 
% The common influence source in the model is 0.5*wi(3,t)
% Verify that the p-value (pValue) is or should be very small, i.e. $p << 0.001$.

%% Concluding remarks 
% * Check the output obtained in this simulation with the one presented in Gourevitch et al.'s article
%
% * Try different sample size to see its influence on estimation and 
%   threshold level. Try for instance, in the [80 to 2000] range to see 
%   how the threshold value in PDC and DC varies.
%   
