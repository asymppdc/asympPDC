%% COMPARE AND TEST FAA X ORIGINAL ASYMPTOTIC ROUTINES
%       This script test and compare the accuracy and speed of FastAsympAlg
%       and the original implementation of asympPDC Package's asymp_pdc and
%       asymp_dtf routines.
%
%% Test model used to generate data sample
%          7-dimension VAR[2] model with loop and feedback borrowed from:
%
%          Baccala & Sameshima (2001b). Overcoming the limitations of correlation
%          analysis for many simultaneously processed neural structures,
%          Progress in Brain Research, 130:33--47.
%
%                <http://dx.doi.org/10.1016/S0079-6123(01)30004-3>

%% Data sample generation

clear; clc; close all;
format compact; format short

% Seeds the random number generator
rng('shuffle')

if ~isOctave()
   v = ver('MATLAB');
   versionNumberMatlab = str2double(v.Version);
else
   versionNumberMatlab = 8.4; % Does Octave have same compatibility as ML R2014b?!
end
% For MATLAB release history see: <https://en.wikipedia.org/wiki/MATLAB>
% Release      Version
%  2014b         8.4
%  2017a         9.2
%  2020a         9.9
%  2021b         9.11

alpha  = 0.001;  % Significance level for hypothesis tests
nFreqs = 128;   % Number of frequency points

%         [ PDC  gPDC  iPDC   DTF   DC   iDTF ]
flgPlot = [  1     0     0     1     0     0  ]; % Choose measures to xplot

if exist('report_faa_x_asymp_PDC.txt') == 2
   delete('report_faa_x_asymp_PDC.txt')
end

diary report_faa_x_asymp_PDC.txt

nDiscard = 50000;     % number of burn in points of simulation
nPoints  = 2000;      % number of sample size

u = fbaccala2001b_model1_feedback( nPoints, nDiscard );
% Model used for simulation
strID = 'Baccala & Sameshima (2001b) Model I - 7-dimensional VAR[2] w loop+feedback';  

%chLabels = []; % or
chLabels = {'x_1';'x_2';'x_3';'x_4';'x_5';'x_6';'x_7'};

fs = 1; % Sampling frequency

%%
% Plotting parameters
flgPrinting = [1 1 1 2 3 0 0]; % overriding default setting
flgScale = 3;
flgMax = 'tci';
flgSignifColor = 3;
flgColor = [1];
w_max = fs/2;

flgPrintResults = 1; % GCT and/or iGCT result printing flag

%%
% Data pre-processing: detrending and normalization options

flgDetrend = 1;     % Detrending the data set flag
flgStandardize = 0; % Standardization flag 

[nChannels,nSegLength] = size(u);

if nChannels > nSegLength
   u = u.';
   [nChannels,nSegLength] = size(u);
end
if flgDetrend
   for i=1:nChannels, u(i,:) = detrend(u(i,:)); end
   disp('Time series were detrended.');
end
if flgStandardize
   for i=1:nChannels, u(i,:) = u(i,:)/std(u(i,:)); end
   disp('Time series were scale-standardized.');
end

disp('Running MVAR estimation routine.')

maxIP = 30;         % maximum model order to consider.
alg = 1;            % 1: Nutall-Strand MVAR estimation algorithm
criterion = 1;      % 1: AIC, Akaike Information Criteria

tic
[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);
tmvar = toc; 

disp(['Number of channels = ' int2str(nChannels) ' with ' ...
    int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%%
% Testing the adequacy of MAR model fitting through Portmanteau test

h = 20; % testing lag
MVARadequacy_signif = 0.05; % VAR model estimation adequacy significance
                            % level
aValueMVAR = 1 - MVARadequacy_signif; % Confidence value for the testing

[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueMVAR,h,...
                                                      flgPrintResults);
%%
% Granger causality test (GCT) and instantaneous GCT

gct_signif  = alpha;  % Granger causality test significance level
igct_signif = alpha;  % Instantaneous GCT significance level
[Tr_gct, pValue_gct] = gct_alg(u,A,pf,gct_signif,flgPrintResults);

%%
% Loop thru DTF or PDC, and their various metrics: euc, diag and info.
%

vmeasure = {'PDC', 'DTF'};
vmetric  = {'euc', 'diag', 'info'};
kflgPlot = 0; kPlot = 0;

for kmeasure = 1:length(vmeasure)  % PDC or DTF loop
   measure = vmeasure{kmeasure};

   for kmetric = 1:length(vmetric) % metrics type loop
      metric = vmetric{kmetric};
      if kPlot, flgPrinting = [1 1 1 2 0 0 0];  end
      kflgPlot = kflgPlot + 1;

%%
% PDC or DTF estimation using original asymp_pdc or asymp_dtf routine
%

      if isempty(regexp(lower(measure),'pdc', 'once'))
         switch metric
            case 'euc'
               measure = 'DTF';
            case 'diag'
               measure = 'gDTF';
            case 'info'
               measure = 'iDTF';
         end
         % disp('asymp_dtf')
         tic
         d = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Estimate DTF and asymptotic statistics
         toriginal = toc;
      else
         switch metric
            case 'euc'
               measure = 'PDC';
            case 'diag'
               measure = 'gPDC';
            case 'info'
               measure = 'iPDC';
         end
         % disp('asymp_pdc')
         tic
         d = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
         toriginal = toc;
      end

      d.Tragct = Tr_gct;
      d.pvaluesgct = pValue_gct;

      if  flgPlot(kflgPlot)
         if isempty(strfind(measure,'DTF')), flgPrinting(5) = 3; end % Print GCT p-values for PDC
         [hxlabel,hylabel] = xplot(strID,d,flgPrinting,fs,w_max,chLabels, ...
                                   flgColor,flgScale,flgMax,flgSignifColor);

         strTitle1 = ['Original Asymp ' measure];
         strTitle2 = [' [ N=' int2str(nSegLength) 'pts; IP=' int2str(d.p)  ' ]'];
         strTitle  = [strTitle1 strTitle2];
         xplot_title(alpha,metric,measure,strTitle);
         kPlot = kPlot + 1;

      end

%%
% DTF or PDC estimation through Fast Asymptotic Algorithm
%

      NNN = 1;
      if versionNumberMatlab > 9.2
%          disp('FastAsympAlg Original')
         tic
         for kkkk=1:NNN
            [cf,df] = FastAsympAlg_orig(u,A,pf,1:nFreqs,measure,[],alpha);
         end
         tfast = toc;
      else
%          disp('FastAsympAlg adapted')
         tic
         for kkkk=1:NNN
            [cf,df] = FastAsympAlg(u,A,pf,1:nFreqs,measure,[],alpha);
         end
         tfast = toc;
      end
        
      df.Tragct = Tr_gct;
      df.pvaluesgct = pValue_gct;

      if flgPlot(kflgPlot)
         if isempty(strfind(measure,'DTF')), flgPrinting(5) = 3; end % Print GCT p-values for PDC
         [hxlabel,hylabel] = xplot(strID,df,flgPrinting,fs,w_max,chLabels, ...
                                   flgColor,flgScale,flgMax,flgSignifColor);

         strTitle1 = ['FastAsympAlg: ' measure];
         strTitle2 = [' [ N=' int2str(nSegLength) 'pts; IP=' int2str(d.p) ' ]'];
         strTitle = [strTitle1 strTitle2];
         xplot_title(alpha,metric,measure,strTitle);

         kPlot = kPlot + 1;

      end

%%
% Comparing estimates between original asymp_PDC package asymp_pdc/asymp_dtf
% function and 'Rezaei et al.(2022)'s FastAsympAlg implementations
%

      deltaFAACIupper = cf.CIupperbound - cf.Phi;
      if isempty(regexp(lower(measure),'pdc', 'once'))
         deltaAsympCIupper = d.ci2 - d.dtf2;
         L = d.dtf2;
      else
         deltaAsympCIupper = d.ci2 - d.pdc2;
         L = d.pdc2;
      end
      tol = 1e-12; % Tolerance for comparison.

      fprintf(['\n\n' repmat('=',1,100) '\n'])
      fprintf(['     COMPARING ' measure ...
                        '2 RELATED ESTIMATES BETWEEN THE ORIGINAL AND FAA ROUTINES \n'])
      fprintf([repmat('=',1,100) '\n'])

      fprintf(1,['\n**  i. MVAR processing time:'])
      fprintf(1,['\n               MVAR = %7.5g s.\n'], tmvar)
      fprintf([repmat('-',1,100) '\n'])
      
      fprintf(1,['\n** ii. Asymptotic routines processing times:'])
      if ~isempty(strfind(measure,'PDC'))
         fprintf(1,['\n           asymp_PDC = '])
      else
         fprintf(1,['\n           asymp_DTF = '])         
      end
      fprintf(2,['%10.5f s;'], toriginal)

      fprintf(1,['\n        FastAsympAlg = '])
      fprintf(2,['%10.5f s;'], tfast)
      
      fprintf(1,['\n       Gain in speed = '])
      fprintf(2,['%7.2f times.\n'], toriginal/tfast)
      fprintf([repmat('-',1,100) '\n'])
      
      fprintf(1,['\n   tol = %4.3g; nFreqs = %1.0f; \x03b1 = ' ...
                                             num2str(alpha) '.\n'], tol, nFreqs)
      fprintf([repmat('-',1,100) '\n'])

%%  Total sum of all absolute differences of PDF2 or DTF2 estimates 
      fprintf(['\n** 1. Comparing ' measure ...
                                  '2 values: Fast x Original ASYMP routines:\n'])

      maxMeasureDiffs = max(max(max(abs(cf.Phi(:,:,:) - L(:,:,:)))));
      
      if maxMeasureDiffs > tol, kwarning = 2; else kwarning = 1; end
      fprintf(kwarning,['      * Maximum absolute difference of ' ...
                                                       measure '2 estimates: '])
      fprintf(2,['          %7.5g. \n'], maxMeasureDiffs)

%       totalMeasureDiffs = sum(sum(sum(abs(cf.Phi(:,:,:) - L(:,:,:)))));
%       if totalMeasureDiffs > tol, kwarning = 2; else kwarning = 1; end
%       fprintf(kwarning,['      * Total sum of all absolute differences of ' ...
%                             measure '2 estimates: %7.5g. \n'], totalMeasureDiffs)

%       assert(maxMeasureDiffs <= tol, [ measure ...
%                                         '2 estimates values not as expected.'])

%%   Total sum of absolute differences of {PDC2 or DTF2} THRESHOLD estimates 
      fprintf(['\n** 2. Comparing threshod values for ' measure ...
                                        '2: Fast x Original ASYMP routines:\n'])

      maxMeasureThresDiffs = max(max(max(abs(cf.Threshold(:,:,:) ...
                                                              - d.th(:,:,:)))));
      if maxMeasureThresDiffs > tol, kwarning = 2; else kwarning = 1; end
      fprintf(kwarning,['      * Maximum absolute difference of ' ...
                                             measure '2 threshold estimates: '])
      fprintf(2,['%7.5g. \n'], maxMeasureThresDiffs)

%       totalMeasureThresDiffs = sum(sum(sum(abs(cf.Threshold(:,:,:) ...
%                                                               - d.th(:,:,:)))));
%       if totalMeasureThresDiffs > tol, kwarning = 2; else kwarning = 1; end
%       fprintf(kwarning,['      * Total sum of absolute differences of ' ...
%                       measure '2 threshold: %7.5g. \n'], totalMeasureThresDiffs)

%       assert(maxMeasureThresDiffs <= tol, ['Threshold ' measure ...
%                                              ' estimates not as expected.'])

%%   Total sum of absolute differences of {PDC2 or DTF2} CI estimates 
      fprintf(['\n** 3. Comparing CI for ' measure ...
                                         '2: Fast x Original ASYMP routines:\n'])
       
      maxMeasureCIDiffs = max(max(max(abs(cf.CIupperbound(:,:,:) ...
                                                             - d.ci2(:,:,:)))));
      if maxMeasureCIDiffs > tol, kwarning=2; else kwarning=1; end
      fprintf(kwarning,['      * Maximum absolute difference of ' ...
                                                   measure '2 CI estimates: '])
      fprintf(2,['       %7.5g. \n'], maxMeasureCIDiffs)

%       totalsumMeasureCIDiffs = sum(sum(sum(abs(cf.CIupperbound(:,:,:) ...
%                                                              - d.ci2(:,:,:)))));
%       if totalsumMeasureCIDiffs > tol, kwarning=2; else kwarning=1; end
%       fprintf(kwarning,['      * Total sum of absolute differences of ' ...
%                             measure '2 CI: %7.5g. \n'], totalsumMeasureCIDiffs)

%       assert(maxMeasureCIDiffs <= tol, ['CI Upper Bound for ' measure ...
%                                                ' estimates not as expected.'])
      
   end
end

fprintf(['\n' repmat('=',1,100) '\n'])

diary off
if ~isOctave() && sum(flgPlot)
   if kPlot < 5
      tilefigs2([],[],2,2)
   else
      tilefigs2;
   end
end

%% Logs:
% 
