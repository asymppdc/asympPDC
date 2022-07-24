%% Andrews and Herzberg Sunspot-Melanoma 1936-1872 Series
%
%% Data File: skin.dat
% This data is from Andrews and Herzberg.
% D. F. Andrews, A. M. Herzberg. (1985) Data: A Collection of Problems from 
% Many Fields for the Student and Research Worker. Springer, New York.
%
% The aetiology of melanoma is complex and may include the influences
% of trauma, heredity and hormonal activity. In particular, exposure
% to solar radiation may be involved in the pathogenesis of melanoma.
% Melanoma is more common in fair-skinned individuals and most frequent
% in skin sites exposed to the sun. In white populations melanoma is
% more common in areas closer to the equator where the intensity of
% solar radiation is higher. Data from various parts of the world suggest
% that the incidence of melanoma is increasing. The data below, giving
% age-adjusted melanoma incidence, are from the Connecticut Tumor Registry
% from 1936-1972. Connecticut has the longest record of state
% population-based cancer statistics in the United States of America.
% The data also includes the sunspot relative number. Houghton, Munster
% and Viola (1978) have shown that the age-adjusted incidence rate for
% malignant melanoma in the state of Connecticut has risen since 1935 and
% that superimposed on the rise are 3-5 year periods in which the rise
% in the rate of incidence is excessive. These periods have a cycle of
% 8-11 years and follow times of maximum sunspot activity. The relationship
% between solar cycles and melanoma supports the hypothesis that melanoma
% is related to sun exposure and provides evidence that solar radiation
% may trigger the development of clinically apparent melanoma. The columns
% are the year, male incidence, total incidence, and sunspot relative
% index. The incidence are rates per 100,000.
%
% Below is the contents of the file skin.dat:
% Year % 1936-1972
% ANNUAL MALE MELANOMA INCIDENCE(AGE-ADJUSTED PER 10**5) CONNECTICUT
% ANNUAL TOTAL MELANOMA INCIDENCE(AGE-ADJUSTED PER 10**5) CONNECTICUT
% ANNUAL SUNSPOT RELATIVE NUMBER

%%============= Causality Analysis using PDC ==============================
% In this example, ANNUAL SUNSPOT RELATIVE NUMBER and TOTAL MELANOMA
% INCIDENCE(AGE-ADJUSTED PER 10**5) in the state of CONNECTICUT will be
% considered

clear; clc; format compact
close all

format short
warning('off'); more off

disp('======================================================================');
disp('     Andrews and Herzberg''s Sunspot and Melanoma 1936-1972 Data');
disp('                 Sunspot --> Melanoma or other way?');
disp('======================================================================');

flgDataPlot = 1;

%% Retrieving data set for analysis
%

x=sunmeladat([4 3]); chLabels={'Sunspot';'Melanoma'};

u = x;

nPoints = length(u);

%%
% Plotting sunspot and melanoma series
  
if flgDataPlot,
   year=sunmeladat(1);  % Year variable
   h1=figure;
   set(h1,'NumberTitle','off','MenuBar','none', ...
      'Name','Andrews & Herzberg''s Sunspot-Melanoma Series (1936-1972)')
   h = subplot(311); plot(year,u(:,1),'o:r'); title('Sunspot series'); grid
   set(h,'XLim', [1935 1975], ...             %'YLim',ylim, ...
      'XTick',[1935:5:1975], ...
      'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ']);

   h = subplot(312); plot(year,u(:,2),'.:'); 
                     xlabel(''); 
                     title('Melanoma series');grid
   set(h,'XLim', [1935 1975], ...             %'YLim',ylim, ...
      'XTick',[1935:5:1975], ...
      'XTickLabel',[' ';' ';' ';' ';' ';' ';' ';' ';' ']);
   h = subplot(313); plot(year,detrend(u(:,2)),'.:'); xlabel('year');
   set(h,'XLim', [1935 1975])
   title('detrended melanoma series')
   grid; hold on; pause(1)

   h2=figure;
   set(h2,'NumberTitle','off','MenuBar','none', ...
      'Name','Andrews & Herzberg''s Sunspot-Melanoma Series (1936-1972)')
   u2=u(:,2); u2=detrend(u2);
   h = subplot(211); plot(year,(u2-mean(u2))/std(u2),'b*:',...
      year,(u(:,1)-mean(u(:,1)))/std(u(:,1)),'ro-');
   set(h,'XLim', [1935 1975])
   title('Melanoma & Sunspot: Standardized series')
   legend('Melanoma','Sunspot');grid

   [uxcorr,lag]=xcorr(u2,u(:,1),'coeff');
   subplot(212); plot(lag,uxcorr,'o--'); xlabel('lag(year)');
   title('Normalized cross-correlation function'); axis([-20 20 -1 1]); 
   grid; hold on
   %if ~is_octave(), tilefigs; end;
%%
% Note cross-correlation peak at 2-year time lag.
end;

%%
% Data pre-processing: detrending and standardization options

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

alpha = 0.05;
fs= 1;
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
[Tr_gct, pValue_gct, Tr_igct, pValue_igct] = gct_alg(u,A,pf,gct_signif, ...
                                              igct_signif,flgPrintResults);
pause(3);

%% Original PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_pdc.m or issue
%
%   >> help asymp_pdc 
%
% command for more detail.
nFreqs = 64;
metric = 'info';  % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
alpha = 0.01;
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics


%%
% PDCn Matrix Layout Plotting
                
% ---------------Plotting options flag setting-----------------------------
%           [1 2 3 4 5 6 7]
flgPrinting=[1 1 1 2 2 0 1];
%            | | | | | | 7 Spectra(0: w/o SS; 1: Linear; 2: log-scale)
%            | | | | | 6 Coherence
%            | | | | 5 Plot lower confidence limit
%            | | | 4 Plot upper confidence limit
%            | | 3 Significant PDC(w) in red line
%            | 2 Patnaik threshold level in black dashed-line
%            1 PDC in green line
%--------------------------------------------------------------------------
flgColor = 0;
w_max=fs/2;
flgScale = 3;
flgMax = 'tci';
flgSignifColor = 1;

for kflgColor = flgColor,
%    h=figure;
%     set(h,'NumberTitle','off','MenuBar','none', ...
%        'Name','Andrews & Herzberg''s Sunspot-Melanoma Series (1936-1972)')
strBarTitle = 'Andrews & Herzberg''s Sunspot-Melanoma Series (1936-1972)';
[hxlabel,hylabel] = xplot(strBarTitle,c,flgPrinting,fs,w_max,chLabels, ...
                                 flgColor,flgScale,flgMax,flgSignifColor);
   xplot_title(alpha,metric,'pdc',strBarTitle);
end;


pause(4)
%%
% One interesting aspect of asymptotic statistics you may notice is that it
% gives proper and consistent inference for all three forms of
% PDC (PDC|gPDC|iPDC) and DTF (DTF|DC|iDTF). 
%
% This is of one indication of the robustness of the algorithm based on
% asymptotic statistics implemented in the AsympPDC Package. 

%% Generalized PDC estimation and analysis parameters
%
% gPDC analysis results are saved in *d* structure.
% See asymp_pdc.m or issue
%
%   >> help asymp_pdc 
%
% command for more detail.
alpha = 0.05;
metric = 'info'; % euc  = original PDC;
                % diag = generalized PDC or gPDC;
                % info = information PDC or iPDC
d = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics


%% 
% PDCn Matrix Layout Plotting
flgColor = 0;
flgScale = 2;
flgMax = 'tci';
flgSignifColor = 3;

% h=figure;
% set(h,'NumberTitle','off','MenuBar','none', ...
%    'Name','Andrews & Herzberg''s Sunspot-Melanoma Series (1936-1972)')
strBarTitle = 'Andrews & Herzberg''s Sunspot-Melanoma Series (1936-1972)';
[hxlabel,hylabel] = xplot(strBarTitle,d,flgPrinting,fs,w_max,chLabels, ...
                                 flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric,'pdc',strBarTitle);

%%
%
% Note also that peak of gPDC in frequency coincides with the peack or
% Sunspot and melanoma series power spectra, which corresponds to 
% approximately 11-year Sun's activity cycle period. 

