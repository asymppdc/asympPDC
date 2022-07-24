%PDC analysis getting started template file
%
% Edit this file to analyze your data. You might want to choose analysis
% parameters followed by comment containing "<***>". Check bellow.
%
% Some important input and output parameters and variables:
% input:
%        u     - data in rows
%        fs    - Sampling frequency
%        maxIP - externally defined maximum IP
%        alg   - for algorithm (1: Nutall-Strand),(2: mlsm) ,
%                              (3: Vieira Morf),  (4: ARfit)
%        criterion - for AR order selection =>
%                                   1: AIC; 2: Hanna-Quinn; 3: Schwarz;
%                                   4: FPE, 5: fixed order in MaxIP
%        alpha  -  PDC test significance level
%
% output: c struct variable from asymp_pdc or asymp_dtf routine.
%
%         c.pdc2 - squared original/generalized/information PDC
%         c.th -  threshold level by Patnaik approximation
%         c.pdc_th - above threshold pdc values otherwise equal NaN
%         c.ic1,c.ic2 - superior and inferior confidence interval
%         c.p - VAR model order
%         c.SS - Power spectrum
%         c.coh - coherece function


%===========================#
% Times series for analysis /
%===========================#
% u     - data in rows.
%         The variable u must contain the time series
%         For instance, if flgExample = 1 the template file will analyze a
%         5-variable Gaussian independent noise processes.

format compact
clear all; clc

% Choose Example 1 == Five independent Gaussian random variables model
%                2 == Sunspot-melanoma time series
%                3 == Baccala & Sameshima (2001) 5-variable linear model
%                4 == Takahashi(2008) Thesis' example model                
flgExample = 4;

disp('======================================================================');
disp('============= PDC analysis getting started template ==================')
disp('======================================================================');

switch flgExample
   case 1,
      u=randn(2000,5);    %<***> Example (1)
      disp('            Random Independent Process with 5 variables')
      disp('======================================================================');

   case 2
      u=sunmeladat([4 3]); %<***> Example (2)
      disp('     Andrews and Herzberg''s Sunspot and Melanoma Example');
      disp('                 Sunspot --> Melanoma or other way?');
      disp('======================================================================');
   case 3
      u=baccala2001a_ex5data(2000);
   case 4,
      u=takahashi_thesis_dat(200);
   otherwise
      error('Wrong example selection.')
end;

fs = 1; %<***> Sampling frequency

[nSegLength,nChannels]=size(u);
if nSegLength < nChannels, error('The data might be transposed.'); end;

%===========================#
%  Channel identification   /
%===========================#

switch flgExample
  case 1,
    chLabels  = {'x_1';'x_2';'x_3';'x_4';'x_5'}; %<***> Example (1)
    strTitle2 = 'Five independent Gaussian noises '; %Title info
  case 2
    chLabels  = {'Sunspot';'Melanoma'}; %<***> Example (2)
    strTitle2 = 'Sunspot-Melanoma 1936-1972 ';
  case 3,
    chLabels  = [];                     %<***> Example (3)
    strTitle2 = 'Five-variable example ';
  case 4,
    chLabels  = {'X'; 'Y'; 'Z'};        % Takahashi thesis example
    strTitle2 = 'Takahashi 2008 (Thesis) example';
end;

flgLabels = ~isempty(chLabels);
if flgLabels,
  if nChannels ~= max(size(chLabels))
    error('Numbers of labels and channels do not match.')
  end;
end;

%==========================================#
%       Pre-processing Action flags        /
%==========================================#
flgDetrend = 1;     %<***> Usually it's recommended to detrend the time series.

flgStandardize = 0; %<***> For PDCn estimation normalization has no effect.
if flgStandardize,
  disp('Be aware that the data standardization does not affect the generalized')
  disp('   and information PDC/DTF estimates nor its statistical results, ')
  disp('  so that data standardization is not necessary in these measures.')
end;

%===========================#
%    Analysis parameters    /
%===========================#
nFreqs = 128; %<***> number of points on frequency scale; 
              %      recommended to use either 64 or 128.

metric = 'info';
%        metric   'euc'  - Euclidean     -> original PDC or DTF;
%                 'diag' - diagonal      -> gPDC or DC;
%                 'info' - information   -> iPDC or DTF;

maxIP = 30; % maxIP - externally defined maximum IP %<***>

%===========================#
%    MAR algorithm          /
%===========================#
% Choose one of algorithm for MAR estimation
% alg   - for algorithm  (1: Nutall-Strand),(2: mlsm) ,
%                        (3: Vieira Morf),  (4: QR artfit)
alg = 1; %<***> Nuttall-Strand (alg=1) algorithm, it seems to be a good 
         %      and robust method.

%============================#
%MAR order selection criteria/
%============================#
% criterion - for AR order choice
%  1: AIC; 2: Hanna-Quinn; 3: Schwarz or BIC;
%  4: FPE, 5: fixed order in MaxIP
criterion = 1; %<***> AIC, Akaike information criterion (Our preferred one)

%==================
alpha = 0.01;        %<***> Significance level for PDC null hypothesis 
                     % testing,
                     %
                     % Note: if alpha == 0, no asymptotic statistics 
                     % calculation is performed and ASYMP_PDC (see bellow) 
                     % will only returns PDC. This option is interesting 
                     % if you want faster PDC calculation.
                     %
gct_signif = alpha;  % Granger causality test significance. Choose other 
                     % value if you have good reason for adopting different 
                     % one from frequency domain statistical testing. 
igct_signif = alpha; % Instantaneous Granger causality test significance level.
VARadequacy_signif = 0.05; % VAR model adequacy significance level

%==========================================================================
%                    Plotting options
%==========================================================================

flgScale = 2; % 1: [0 1] / {if max(PDC2/DTF2) > 1}:[0 max(PDC2/DTF2)]
%               2: [0 {if max(PDC/DTF2) > 1}:max(PDC/DTF2)]/[0 1]/[0 .1]/[0 .01]
%                                         based on flgMax (PDC2/DTF2/Thr/CI/all)
%               3: [0 max(PDC2/DTF2/Thr/CI/all)]
%                                         based on flgMax (PDC2/DTF2/Thr/CI/all)
%               4: [0 {max(PDC2/DTF2/Thr/CI/all) or round to {0.01 0.1 1.0}]
%
flgMax = 'all'; % {'PDC'|'DTF'|'Thr'|'CI'|'TCI'|'all'} measure used as upper limit
%                for y-axis scales:
%                    PDC|DTF - PDC2/DTF2 maximum value;
%                        Thr - Patnaik threshold maximum value;
%                        CI  - maximum upper confidence interval value;
%                        TCI - threshold or CI max value;
%                        all - plotted either max([PDC2/DTF2, Thr, CI]) value;
%               See also flgScale.

flgSignifColor = 3; % 0: black line
                    % 1: black / gray -> significant /not signif PDC/DTF
                    % 2: red  /  gray ->      "
                    % 3: red  / green ->      "
                    % 4: red  / black ->      "
                    % 5: black / green

flgColor = [0 1]; % Plotting option for automatic scaling for small PDC
                  % values.
                  % if flgColor = 0, y-axis scale = [0 1]
                  % elseif flgColor = 1, the pdc_xplot routine rescales 
                  % the y-axis automatically according to the following 
                  % rules:
                  %   if .01<=PDC(f) < .1 background-color = light-blue,
                  %                          so that y-axis scale = [0 .1]
                  %   elseif PDC(f) < .01 background-color = light-purple
                  %                          and y-axis = [0 .01];
                  % for flgColor=[0 1], both lay-outs are plotted.


%                 1 2 3 4 5 6 7
flgPrinting   =  [1 1 1 2 2 0 1]; % Example: Plot everything, except coh2.
%           blue  | | | | | | 7-- {0:5} Spectra (0: w/o; 1: Linear; 2: Log;
%                 | | | | | |           3: PDC2; 4: Linear normalized;
%                 | | | | | |           5: Log spectra + PDC2)
%           gray  | | | | | 6-- {0:1} Coh2 (0: w/o Coh2; 1: w Coh2)
%    dark-purple  | | | | 5-- {0:3} Print GCT p-values and dot-mark significant
%  or dark-green  | | | |           connectivity channel-pair (0: w/o;
%                 | | | |           1: p-values; 2: dot-mark +GCT;
%                 | | | |           3: p-values + dot-mark significant GCT)
%    dashed-blue  | | | 4-- {0:4} Confidence interval (0:w/o; 1: Dashed-lines;
%                 | | |           2: Shaded-plot; 3: Error-bar 1; 4: Error-bar 2
%            red  | | 3-- {0:1} Significant PDC2|DTF2 in red lines
%   dashed-black  | 2-- {0:1} Patnaik threshold level in black dashed-lines
%          green  1-- {0:1} PDC2/DTF2 in green lines or black w/o statistics,
%                           see flgSignifColor bellow for line color selection.


w = fs*(0:(nFreqs-1))/2/nFreqs;
w_max = fs/2; %<***> Usually half of sampling frequency = Nyquist frequency


%==========================================================================
%==========================================================================
%        WARNING: BELOW THIS LINE PROBABLY YOU MIGHT NOT WANT TO EDIT,
%            UNLESS YOU NEED TO CUSTOMIZE YOUR ANALYSIS ROUTINE.
%==========================================================================

%==========================================================================
%                    Detrend and standardization options
%==========================================================================


% Determine time series length.
[nChannels,nSegLength]=size(u);
if nChannels > nSegLength, u=u.'; 
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



%==========================================================================
% Additional info for title (optional)

strTitle1 = [];

%==================
switch alg
  case 1
    disp('VAR estimation using Nutall-Strand algorithm.')
  case 2
    disp('VAR estimation using least-squares estimator.')
  case 3
    disp('VAR estimation using Vieira-Morf algorithm.')
  case 4
    disp('VAR estimation using QR-Arfit algorithm.')
end;

%============================#
%MAR order selection criteria/
%============================#
switch criterion
   case 1
      disp('Model order selection criteria: AIC.')
   case 2
      disp('Model order selection criteria: Hanna-Quinn.')
   case 3
      disp('Model order selection criteria: Schwarz (BIC).')
   case 4
      disp('Model order selection criteria: FPE.')
   case 5
      disp('Model order selection criteria: fixed order in maxIP.')
   otherwise
      error('Model order selection criteria: NOT IMPLEMENTED YET.')
end;

%==========================================================================
%                            VAR model estimation
%==========================================================================
[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);


disp(['Number of channels = ' int2str(nChannels) ' with ' ...
  int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%==========================================================================
%    Testing for adequacy of MAR model fitting through Portmanteau test
%==========================================================================
   h = 20; % testing lag
   aValueVAR = 1 - VARadequacy_signif;
   flgPrintResults = 1;
[Pass,Portmanteau,st,ths]=mvarresidue(ef,nSegLength,IP,aValueVAR,h,...
                                                          flgPrintResults);

%==========================================================================
%         Granger causality test (GCT) and instantaneous GCT
%==========================================================================
   flgPrintResults = 1;
[Tr_gct, pValue_gct, Tr_igct, pValue_igct] = gct_alg(u,A,pf,gct_signif, ...
                                                         flgPrintResults);

%==========================================================================
%            PDC, threshold and confidence interval calculation.
%==========================================================================

% if alpha == 0, no asymptotic statistics is performed. ASYMP_PDC returns
% only the PDC. This option is much faster!!
tic
  c=asymp_pdc(u,A,pf,nFreqs,metric,alpha);
toc
%                    or 
% c=asymp_dtf(u,A,pf,nFreqs,metric,alpha);

%  Two lines bellow are unnecessary as asymp_dtf and asymp_pdc already 
%  calculate SS and coh.

%Adding further analysis details in the figure title.
%strTitle3 = ['[N=' int2str(nSegLength) '; IP=' int2str(c.p) ']'];
% or

% strTitle3 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p) '; ' ...
%    datestr(now) ']'];
strTitle3 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p) ']'];

% or leave emptied: strTitle3=[];

%==========================================================================
%              Matrix Layout Plotting of the Analysis Results
%==========================================================================

w_max = fs/2;
strTitle = [strTitle1 strTitle2 strTitle3];
strWindowName = 'PDC/DTF/GCT Analysis Template Example';

% The following "for loop" through flgColor values, 0 and 1, and yields a
% pair of plots, one without and other with color scale rearrangement option.
% Value range of PDC and Coherence is from [0 1], but sometimes the maximum 
% peak value is small (<0.1), or even smaller, (<.01), so in these cases it
% might be interesting to have a plot with finer smaller y-axis scale. The
% white-background plot indicates full-scale [0 1] y-axis, while
% light-blue-background stands for intermediate [0 .1] scaling and
% light-purple-background shows very fine detail of small, usually not
% significant PDCs. Try flgColor = 0 or 1, or both [0 1].

for kflgColor = flgColor,
   [h1,~,~] = xplot(strWindowName,c,...
                    flgPrinting,fs,w_max,chLabels,kflgColor);
  
   [h2,~,~] = xplot(strWindowName,c,flgPrinting,fs,w_max,chLabels, ...
                                 kflgColor,flgScale,flgMax,flgSignifColor);
   
% The title suplabel command should (not sure) follow the dtf_xplot routine
% In MacOS X, for flgPrinting(7) = 4 or 5, the main diagonal plotting
% gets misaligned if suplabel with 't' option is used more than once.
   xplot_title(alpha,metric,'pdc',strTitle);
end;


%======================= xplot ========================================
%Plot legend:  Blue lines on the main diagonal = Power spectra;
%              Black dashed lines are Patnaik threshold for dtfn;
%              Green lines = non significant dtfn;
%              Red lines = significant dtfn;
%              Light-gray lines = coherence function.
%
% Notes:       a.The main diagonal of matrix layout contains power spectra.
%              b.Coherences are symmetric, e.g., 
%                   Coh_{Sunspot,Melanoma}(f) = Coh_{Melanoma,Sunspot}(f).
%              c.dtfn is asymmetric relation, and the dtfn graphics should
%              be read as if the flow of information is been from the
%              x-axis variable toward y-axis variable.
%
%              For sunspot and melanoma example, one only sees significant
%              dtfn from Sunspot to Melanoma, which could eventually be
%              interpreted that "Sunspot", or the Sun's activity
%              modulates the incidence of melanoma.
%======================= dtf_xplot ========================================
disp('======================================================================');
disp('===========ANALYSIS_TEMPLATE SUCCESSFULLY FINISHED ===============')
disp('======================================================================');
