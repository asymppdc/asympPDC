%% ANDREWS & HERZBERG (1985)a - Sunspot-Melanoma 1936-1972 Series PDC Demo
%
%% Data File: skin.dat
% This data is from:
%    D. F. Andrews, A. M. Herzberg. (1985) Data: A Collection of Problems from
%    Many Fields for the Student and Research Worker. Springer, New York.
%
% "The aetiology of melanoma is complex and may include the influences of
% trauma, heredity and hormonal activity. In particular, exposure to solar
% radiation may be involved in the pathogenesis of melanoma. Melanoma is more
% common in fair-skinned individuals and most frequent in skin sites exposed to
% the sun. In white populations melanoma is more common in areas closer to the
% equator where the intensity of solar radiation is higher. Data from various
% parts of the world suggest that the incidence of melanoma is increasing. The
% data below, giving age-adjusted melanoma incidence, are from the Connecticut
% Tumor Registry from 1936-1972. Connecticut has the longest record of state
% population-based cancer statistics in the United States of America. The data
% also includes the sunspot relative number. Houghton, Munster and Viola (1978)
% have shown that the age-adjusted incidence rate for malignant melanoma in the
% state of Connecticut has risen since 1935 and that superimposed on the rise
% are 3-5 year periods in which the rise in the rate of incidence is excessive.
% These periods have a cycle of 8-11 years and follow times of maximum sunspot
% activity. The relationship between solar cycles and melanoma supports the
% hypothesis that melanoma is related to sun exposure and provides evidence that
% solar radiation may trigger the development of clinically apparent melanoma.
% The columns are the year, male incidence, total incidence, and sunspot
% relative index. The incidence are rates per 100,000." (Andrews & Herzberg,
% 1985)
%
% Below is the contents of the file skin.dat:
%%
% * Year % 1936-1972
% * ANNUAL MALE MELANOMA INCIDENCE(AGE-ADJUSTED PER 10**5) CONNECTICUT
% * ANNUAL TOTAL MELANOMA INCIDENCE(AGE-ADJUSTED PER 10**5) CONNECTICUT
% * ANNUAL SUNSPOT RELATIVE NUMBER
%
%% Causality Analysis using Partial Directed Coherence (PDC) 
% In this example, ANNUAL SUNSPOT RELATIVE NUMBER and TOTAL MELANOMA
% INCIDENCE(AGE-ADJUSTED PER 10**5) in the state of CONNECTICUT will be
% considered
%
%% Additional References
%
% [0] Andrews, D.F., Herzberg, A.M. (1985). Incidence of Malignant Melanoma
%     After Peaks of Sunspot Activity. In: Data. Springer Series in Statistics.
%     Springer, New York, NY. <https://doi.org/10.1007/978-1-4612-5098-2_33>
% 
% [1] M.V. Viola, A. Houghton, E.W. Munster. Solar cycles and malignant melanoma.
%     Medical Hypotheses 5:153--160, 1979.
%       <https://doi.org/10.1016/0306-9877(79)90067-7>
%
% [2] A. Houghton, E.W. Munster and M.V. Viola. Increased incidence of
%     malignant melanoma after peaks of sunspot activity. The Lancet,
%     311:759--760, 1978.
%       <https://doi.org/10.1016/S0140-6736(78)90869-3>

% (C) Koichi Sameshima & Luiz A. Baccala, 2021. See file license.txt in
% installation directory for licensing terms.

clear; clc; format compact

format short
warning('off'); more off

disp(repmat('=',1,100))
disp('     Andrews and Herzberg''s Sunspot and Melanoma 1936-1972 Data');
disp('                 Sunspot --> Melanoma or other way?');
disp(repmat('=',1,100))

flgDataPlot = 1;

%% Retrieving data set for analysis from sunmeladat.m m-file
%

x=sunmeladat([4 3]); chLabels={'Sunspot';'Melanoma'};
u = x;
nPoints = length(u);

%% Plotting sunspot and melanoma series
%  
if flgDataPlot
   year=sunmeladat(1);  % Year variable
   h1=figure;
   set(h1,'NumberTitle','off','MenuBar','none', ...
      'Name','Andrews & Herzberg''s Sunspot-Melanoma Series (1936-1972)')
   
   h11 = subplot(311); plot(year,u(:,1),'o:r'); title('Sunspot series'); grid
   set(h11,'XLim', [1935 1975], 'XTick',[1935:5:1975], ...
                            'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ']);

   h12 = subplot(312); plot(year,u(:,2),'.:'); 
                     xlabel(''); 
                     title('Melanoma series');grid
   set(h12,'XLim', [1935 1975],'XTick',[1935:5:1975], ...
                            'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ']);
   
   h13 = subplot(313); plot(year,detrend(u(:,2)),'.:'); xlabel('year');
   set(h13,'XLim', [1935 1975])
   title('Detrended Melanoma series')
   grid; hold on;
   shg; pause(3)
   h2=figure;
   set(h2,'NumberTitle','off','MenuBar','none', ...
      'Name','Andrews & Herzberg''s Sunspot-Melanoma Series (1936-1972)')
   u2=u(:,2); u2=detrend(u2);

   h21 = subplot(211); plot(year,(u2-mean(u2))/std(u2),'b*:',...
      year,(u(:,1)-mean(u(:,1)))/std(u(:,1)),'ro:');
   xlabel('year');
   set(h21,'XLim', [1935 1975])
   title('Melanoma & Sunspot: Standardized series')
   legend('Melanoma','Sunspot');grid

   [uxcorr,lag] = xcorr(u2,u(:,1),'coeff');
   subplot(212); h22 = plot(lag,uxcorr,'o'); xlabel('lag(year)');
   title('Normalized cross-correlation function'); axis([-20 20 -1 1]);
   grid; hold on
   plot(2,uxcorr(39),'r.','MarkerSize',15) % Peak of cross-correlation.
   shg; pause(3)
end

%%
% Note the cross-correlation peak at 2-year time lag (filled circle).

%%
% Although for plotting purpose the Sunspot series was standardized, in the
% analyses that follow *the time series were not standardized*, and consequently
% the variance of Sunspot series is much larger than that of the Melanoma's
% series. This may give misleading high *Original PDC* values, although not
% significant, from *Melanoma to Sunspot*, as seen in the following analysis.


%% 
% Data pre-processing: detrending and standardization options


flgDetrend = 1;     % Detrending the data set
flgStandardize = 0; % No standardization
[nChannels,nSegLength] =size(u);
if nChannels > nSegLength, 
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

alpha_c = 0.01;
fs= 1; % Sampling frequency: 1 sample/year.

%% MVAR model estimation
%

maxIP = 4;         % maximum model order to consider.
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
disp(['Number of channels = ' int2str(nChannels) ' with ' ...
    int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%%
% Testing for adequacy of MAR model fitting through the Portmanteau test
h = 20; % testing lag
MVARadequacy_signif = 0.05; % VAR model estimation adequacy significance level
aValueMVAR = 1 - MVARadequacy_signif; % Confidence value for the testing
flgPrintResults = 1;                  % Print results in the Command Window
[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueMVAR,h,...
                                           flgPrintResults);

%% Granger causality test (GCT) and instantaneous GCT
%

gct_signif  = 0.01;  % Granger causality test significance level
igct_signif = 0.01;  % Instantaneous GCT significance level
flgPrintResults = 1; % Flag to print gct_alg.m results on Command Window.
[Tr_gct, pValue_gct,]  = gct_alg(u,A,pf,gct_signif,flgPrintResults);
[Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgPrintResults);

%% Original definition of PDC estimation
%
% PDC analysis results are saved in *c* struct variable.
% See *asymp_pdc.m* for more detail.

nFreqs = 128;
metric = 'euc';  % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
alpha_c = 0.01;
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha_c); % Estimate PDC and asymptotic statistics
c.Tragct = Tr_gct; 
c.pvaluesgct = pValue_gct;

%% Plotting PDC original in Matrix Layout without asymptotic statistics
%

% 
% ---------------Plotting options flag setting-----------------------------
%             [1 2 3 4 5 6 7]
flgPrinting = [1 0 0 0 3 0 1];
%              | | | | | | 7 Spectra(0: w/o SS; 1: Linear; 2: log-scale)
%              | | | | | 6 Coherence
%              | | | | 5 Print GCT p-values
%              | | | 4 Plot upper confidence limit
%              | | 3 Significant PDC(w) in red line
%              | 2 Patnaik threshold level in black dashed-line
%              1 PDC in green line
%--------------------------------------------------------------------------
flgColor = 0;
w_max=fs/2;
flgSignifColor = 0;
flgMax = 'PDC';
flgScale = 2;

strBarTitle = 'Andrews & Herzberg''s Sunspot-Melanoma Series (1936-1972)';
strTitle = 'Sunspot-Melanoma (1936-1972)';

[h1,~, ~] = xplot(strBarTitle,c,...
                          flgPrinting,fs,w_max,chLabels,flgColor, ...
                          flgScale,flgMax,flgSignifColor);
xplot_title(alpha_c,metric,'pdc',strTitle);
shg; pause(3)

%% 
% Note that, in this case, the original PDC from Melanoma (smaller variance) to
% Sunspot (larger variance) estimate is close to unit (large), but not
% significant, while in the opposite direction it is very small and significant.
% The p-values of Granger causality test result are printed above each PDC
% subplot in dark-green (not significant from Melanoma to Sunspot), or in
% dark-purple (significant from Sunspot to Melanoma series, together with a dot
% inside the subplot with the same color).
%
% Proper "scaling" of original PDC estimates could be attained by standardizing
% the time series, i.e scaling to have unit variance and zero mean, or, better,
% using either generalized PDC or information PDC.
%

%% Showing more detail by changing the scale and adding color information
% The following figure shows in more detail the rescaled y-axis of PDC from
% *Sunspot* to *Melanoma*, which one would expect some causal influence of solar
% radiation intensity increase on skin cancer incidence.

flgScale = 2;
flgMax = 'all';
flgSignifColor = 3;

[h2,hxlabel,hylabel] = xplot(strBarTitle,c,flgPrinting,fs,w_max,chLabels, ...
                          flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha_c,metric,'pdc',strTitle);
pause(3)

%%
% One interesting aspect of asymptotic statistics you may notice is that it
% gives proper and consistent inference for all three forms of PDC
% (PDC|gPDC|iPDC) and even for DTF (DTF|DC|iDTF).
%
% The *robustness* is one of such interesting aspects of the analyses based on
% asymptotic statistics implemented in the AsympPDC Package.

%% Generalized PDC estimation and analysis parameters
%
% The gPDC analysis results are saved in *d* struct variable, attained using
% metric = 'diag' (diagonal). See asymp_pdc.m help.

 % Estimate generalized PDC and its asymptotic statistics
metric = 'diag';
alpha_d = alpha_c;
d = asymp_pdc(u,A,pf,nFreqs,metric,alpha_d);

% Assign GCT results to *d* struct variable
d.Tragct = Tr_gct;         
d.pvaluesgct = pValue_gct;

%% Generalized PDC Matrix Layout Plotting
%

flgColor = 0;
[h3,hxlabel,hylabel] = xplot(strBarTitle,d,...
                          flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha_d,metric,'pdc',strTitle);
pause(3)

%%
% The *generalized PDC* provides a better idea  of real _interaction_ between
% *Sunspot* number (an indirect measure of Sun's activity level) and incidence
% of *Melanoma* in the State of Connecticut.
%
% Note that the peak of gPDC from Sunspot to Melanona in frequency coincides
% with the peak of Sunspot and Melanoma series power spectra, which corresponds
% to the 11-year Sun's activity cycle period.


%% Repeating the analysis with Information PDC
%
% iPDC analysis results are saved in *e* struct variable.  
%

% Estimate iPDC and asymptotic statistics
metric = 'info';
alpha_e = 0.01;
e = asymp_pdc(u,A,pf,nFreqs,metric,alpha_e); 
e.Tragct = Tr_gct;
e.pvaluesgct = pValue_gct;

%% iPDC with alpha = 1% analysis Matrix Layout Plotting 
%

flgColor = 0;
flgPrinting=[1 1 1 3 3 1 1];
flgScale = 2;
flgMax = 'all';
flgSignifColor = 3;

[h4,hxlabel,hylabel] = xplot(strBarTitle,e,flgPrinting,fs,w_max,chLabels, ...
                                       flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha_e,metric,'pdc',strTitle);
pause(3)

%%
% Looking closely, the iPDC and gPDC's magnitude patterns, although similar, are
% not equal.

%% Plotting the PDC p-values in the frequency domain 
flgPrinting=[1 1 1 2 3 1 1];
flgScale = 1;
strBarTitle = 'Sunspot-Melanoma Series (1936-1972)'; % Trimmed to reveal p-values

[h5,hxlabel,hylabel] = xplot_pvalues(strBarTitle, c, ...
                                  flgPrinting,fs,w_max,chLabels,flgColor,flgScale);
xplot_title(alpha_c,metric,'p-value iPDC',strTitle);

%%
% In this case $\alpha = 1$% was used, and note the threshold dashed-lines
% located at -2 ($=log_{10}\,0.01$). Only the iPDC from *Sunspot* to
% *Melanoma* is significant in the lower frequency range, which overlaps in
% the frequency scale with peaks of power spectra.

%%
% End of Sunspot-Melanoma series analysis example.

