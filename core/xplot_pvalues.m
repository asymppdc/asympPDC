%% XPLOT_PVALUES
%      Plot pvalues of PDC or DTF connectivity in matrix layout with optional
%      power spectra in the main diagonal.
%
%% Syntax:
%      [hfigure,hxlabel,hylabel] = XPLOT_PVALUES(vBarTitle,c,flgPrinting,
%                                          fs,w_max,chLabels,flgColors,flgScale)
%
%% Input Arguments: 
%
%   strBarTitle:  Title for the figure window bar 
%
%   c.{Tragct,pvaluesgct,SS,Coh,L,Lpatnaik,LTra,
%                                     L2vinf,L2vsup,metric}- results structure
%
%   flgPrinting: [1 1 1 1 1 0 1];
%           blue  | | | | | | 7-- {0:2} Spectra (0: wo; 1: Linear; 2: Log) 
%                 | | | | | 6-- {} Not used
%   dark-purple   | | | | 5-- {0:3} Mark significant GCT pairs + print p-values   
%  or dark-green | | | |          (0: w/o; 1: print p-values; 2: mark +GCT;
%                 | | | |           3: print p-values + mark significant GCT)
%    dashed-blue  | | | 4-- {0:2} Scale of p-values plots 1:linear; 2:log10
%                 | | 3-- {} Not used
%   dashed-black  | 2-- {0:1} Patnaik threshold level in black dashed-line
%                 1-- {} Not used
%
%    fs:        sampling frequency
%    w_max:     frequency scale upper-limit
%    chLabels:  channel identification labels
%    flgColor:  [NOT USED]
%    flgScale:  y-axis scale for p-values: 1: linear [0 1] and log [0 -15];
%                                          2: log y-axis according to 
%                                             max(abs(pvalues_log10)) value
%                                             [0 -5], [0 -10] or [0 -15]
%
%% Output
%     Pretty matrix-layout plots of p-values of |PDC|^2 | |DTF|^2
%     hfigure: figure handle
%     hxlabel,hylabel: label's handles
%
%% Examples: 
%   Calculate c struct results using alg_pdc/asymp_dtf and gct_alg functions.
%
%           xplot_pvalues([],c); % Defaults flgPrinting, fs, w_max and flgColor
%
%           xplot_pvalues([],c,[0 0 0 2 3 0 1],200,100,[],0);
%                        % PDC, threshold, power spectra, coherence
%                        % plotted. fs=200 Hz; default channel label;
%                        % flgColor=0 => no color or rescaling is used.
%
%% Other routines
%  See also  xplot_title, xplot, asymp_pdc, asymp_dtf, FastAsympAlg0, xstats
%                   <xplot_pvalues.html |xplot_pvalues|> |

% (C) Koichi Sameshima & Luiz A. Baccala', 2022. See file license.txt in
% installation directory for licensing terms.

function [hfigure,hxlabel,hylabel] = xplot_pvalues(vBarTitle, c,...
                                flgPrinting,fs,w_max,chLabels,flgColor,flgScale)
knargin = 8;     % Number of input arguments

if nargin < knargin || isempty(flgScale), flgScale = 1; end 
% 1: [0 1] / {if max(PDC2/DTF2) > 1}:[0 max(PDC2/DTF2)]

kylim  = [-0.05  1.075]; % y-axis scale values.

hfigure=figure;
set(hfigure,'PaperOrientation','landscape', ...
            'renderer','painters'); % To get vector figures

if isempty(vBarTitle)
   set(hfigure,'NumberTitle','off','MenuBar','none')
else
   set(hfigure,'NumberTitle','off','MenuBar','none','Name', vBarTitle)
end

if isfield(c,'dtf2')
   L = c.dtf2;       % DTF^2 (N x N x freq)
   if isfield(c,'dtf2_th')
      LTra   = c.dtf2_th; % Significant |DTF|^2 on freq range, otherwise = NaN.
   end
   flgType = 'DTF';
elseif isfield(c,'pdc2')
   L = c.pdc2;       % |PDC|^2 (N x N x freq)
   if isfield(c,'pdc2_th')
      LTra   = c.pdc2_th; % Significant PDC^2 on freq range, otherwise = NaN.
   end
   flgType = 'PDC';
else
   error('Variable c does not hold PDC/DTF analysis results.')
end

if isfield(c,'Tragct')
   Tragct = c.Tragct;
end
if isfield(c,'pvaluesgct')
   pvaluesgct = c.pvaluesgct;
end

if (flgPrinting(1) == 1) && (flgPrinting(5) == 1) && ~isfield(c,'Tragct')
   fprintf(2,['\n*** Warning: to print GCT p-values, assign ''gct_alg()''' ...
              'results to c structure. ***\n'])
   fprintf(2,['\n    Tr_gct, pValue_gct, Tr_igct, pValue_igct] =' ...
                               'gct_alg(u,A,pf,gct_signif,flgPrintResults);'])
   fprintf(2,'\n   then before calling xplot() function execute:')
   fprintf(2,['\n    c.Tragct = Tr_gct;'])
   fprintf(2,['\n    c.pvaluesgct = pValue_gct;\n\n'])
end


SS = c.SS;       % Spectra
Coh2 = c.coh2;     % Coh^2
Lpatnaik = c.th; % Patnaik threshold values for alpha
pvalues = c.pvalues;
metric = c.metric; % Valid optons: "euc", "diag" or "info"
alpha = c.alpha;

[nChannels,~,nFreqs] = size(L); % nChannels:Number of channels/time series;
%                               % nFreqs : number of points on frequency scale
nodesett = 1:nChannels;

% figure scaling parameter adjustment
vChannels = max([10, nChannels]);

%===============================================================================
% Parameters for p-values plotting.
%
N1 = 1; N30 = 30; V1 = 25; V30 = 10;
vMarkerSize = V1 + (vChannels - N1)*(V30 - V1)/(N30 - N1);

vMarkerColor = [ 1       0         1]; % Purple
pValuesNotSignifFontColor = [0 .4 0];
pValuesSignifFontColor = 0.6*vMarkerColor;

N1 = 10; N30 = 24; V1 = 8; V30 = 5; % was 10,5
pValuesFontSize = V1 + (vChannels - N1)*(V30 - V1)/(N30 - N1);

%===============================================================================
% Adjust ticklabels offset gap from axes
%
% (See
% https://www.mathworks.com/matlabcentral/answers/2318-set-position-of-tick-labels
% )
% a=gca;
% a.XRuler.TickLabelGapOffset = -8;
% % negative numbers shift ticklabels downward (positive upward)
% a.YRuler.TickLabelGapOffset = -8;    
% % negative numbers shift ticklabels toward right (negative -> left)

if isOctave()
   f0 = 10; N0 = 4;   % FontSize = 10  for nChannels = 4
   f1 = 6.5; N1 = 20; % FontSize = 6.5 for nChannels = 20
   vTickLabelGapOffset = 0; % YTickLabel distance from y-axis
else
   switch computer %                               ISPC ISUNIX ISMAC ARCHSTR
   %     64-Bit Platforms
   %       PCWIN64  - Microsoft Windows on x64       1     0     0   win64
   %       GLNXA64  - Linux on x86_64                0     1     0   glnxa64
   %       MACI64   - Apple Mac OS X on x86_64       0     1     1   maci64
      case 'PCWIN64'
         f0 = 10; N0 = 4;   % FontSize = 10   for nChannels = 4
         f1 = 6.5; N1 = 20; % FontSize =  6.5 for nChannels = 20
         vTickLabelGapOffset = -1.5; % YTickLabel distance from y-axis  -1.5
      case 'GLNXA64'
         f0 = 9; N0 = 4;    % FontSize = 9   for nChannels = 4
         f1 = 5.5; N1 = 20; % FontSize = 6.5 for nChannels = 20
         vTickLabelGapOffset = -1.0; % YTickLabel distance from y-axis -2.5
      case 'MACI64'
         f0 = 10; N0 = 4;   % FontSize = 10  for nChannels = 4
         f1 = 8;  N1 = 20;  % FontSize =  8  for nChannels = 20
         vTickLabelGapOffset = -0.5; % YTickLabel distance from y-axis  -0.5
      otherwise
         error('computer function not working properly.')
   end
end

%===============================================================================
% FontSize and linewidth scaling parameters adjustment

vFontSize = (vChannels - N0)*(f1 - f0)/(N1 - N0) + f0; % TickLabel Font Size

if vChannels >15
   vLineWidth = 1.25;
   vPatnaikLineWidth = 0.75;
else
   vLineWidth = 1.75;
   vPatnaikLineWidth = 1.25;
end

%===============================================================================
% Set figures size for subplot layout

% Screen dimension in pixel.
set(0,'units','pixels');
sz = get(0,'ScreenSize');

skreenfactor = 2.4; % Ad hoc factor that should 

if ~isOctave()
   % This code is more accurate to get the correct screen size (REF <URL>)
   ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment;
   gd = ge.getDefaultScreenDevice;
   screensize = [gd.getDisplayMode.getWidth gd.getDisplayMode.getHeight];
   
   sz(3) = screensize(1); % width in pixels
   sz(4) = screensize(2); % height in pixels
else
   sz = get(0,'ScreenSize');
end

% Ad hoc check for the presence of multiple monitors in Octave.
khmon = ceil(sz(3)/sz(4)/skreenfactor); % Guessed # of horizontally tiled screens
if khmon < 1, khmon = 1; end

kvmon = round(sz(4)/1000); % Guessed # of stacked screens
if kvmon < 1, kvmon = 1; end
screenratio = sz(3)/sz(4)/khmon;

% Scale figure size on screen according to the monitor resolution
% This has been implemented particularly for Octave.
% Reference monitor has width=sz(3)=1920 pxls
pxwidthscreen = sz(3)/khmon; pxheightscreen = sz(4)/kvmon;

% What follow is a kludge solution to determine figure size in normalized units
% that might allow handling the cases of multiple monitors set up in Octave.
% rwidth  = pxwidth/sz(3)/khmon/kvmon;
% rheight = pxheight/sz(4)/khmon/kvmon;

A4width = 29.7; A4height = 21.0; % A4 landscape paper
A4ratio = sqrt(2); 

% Real figure size for printing or publication in centimeters
width  = min([A4width*vChannels/25  A4width]); 
height = min([A4height*vChannels/25 A4height]); 

%% 
% Determine max figure size in pixels that fit the monitor, which will
% correspond to the A4 max figure size;
if screenratio >= A4ratio % Wider screen than A4 proportion
   pxwidth_max = round(A4ratio * pxheightscreen);
   pxheight_max = pxheightscreen;
else
   pxwidth_max = pxwidthscreen;
   pxheight_max = round(pxwidthscreen/A4ratio);   
end

% Relative figure size on screen scale keeping A4 paper proportion .
rwidth = width/A4width * pxwidth_max/pxwidthscreen/khmon;
rheight = height/A4height * pxheight_max/pxheightscreen;

set(hfigure,'units','centimeters','position',[0 0 width height]);

switch computer
   case 'PCWIN64'
      set(hfigure,'units','normalized', ...
         'position',[(1/khmon-rwidth)/2 (1-1.087*rheight) rwidth rheight]);

   case {'GLNXA64','x86_64-pc-linux-gnu'}
      set(hfigure,'units','normalized', ...
         'position',[(1/khmon-rwidth)/2 (1-1.087*rheight) rwidth rheight]);

   case 'MACI64'
      set(hfigure,'units','normalized', ...
         'position',[(1/khmon-rwidth)/2 (1-1.087*rheight) rwidth rheight]);

   otherwise
      error('computer function not working properly.')
end


%===============================================================================
%         
if nargin <  (knargin-4),  flgPrinting = [0 1 0 0 1 0 1]; end
if nargin <= (knargin-4),  fs = 1; end
if nargin < (knargin-2),   w_max = fs/2; end
if nargin < knargin,       flgColor = 0; end

if w_max > fs/2 + eps
    error(['The parameter w_max should be <= Nyquist frequency,' ...
        'i.e, w_max <= fs/2.'])
end

%flgPrinting = [1 1 1 0 0 0 0];

w = 0:fs/(2*nFreqs):w_max-fs/(2*nFreqs);
nPlotPoints = length(w);
w_min = w(1);

if fs == 1
   str_w_max = '.5';
elseif fs < 1
   str_w_max = sprintf('%0.5g', w_max);
elseif fs < 20
   str_w_max = sprintf('%1.2f', w_max);
   
elseif fs < 200
   if mod(fs,2)
      str_w_max = sprintf('%3.1f', w_max);
   else
      str_w_max = sprintf('%2.0f', w_max);
   end
elseif fs <= 10000
   str_w_max = num2str(w_max); %sprintf('%5.0f', w_max);
else
   str_w_max = 'fs/2';
end

if nargin < (knargin-1) || isempty(chLabels)
    if isfield(c,'chLabels')
        chLabels = c.chLabels;
    else
        chLabels = [];
    end
    if ~isempty(chLabels) && max(size(chLabels)) < nChannels
        error('1 NOT ENOUGH CHANNEL LABELS.');
    end
elseif max(size(chLabels)) < nChannels
    if isfield(c,'chLabels')
        chLabels = c.chLabels
    else
        chLabels = [];
    end
    if ~isempty(chLabels) && max(size(chLabels)) < nChannels
        error('2 NOT ENOUGH CHANNEL LABELS 2.');
    else
        disp('3 NOT ENOUGH CHANNEL LABELS. Default labels assumed.');
    end
end

hxlabel = 0; % x-axis labels' handles
hylabel = 0; % y-axis labels' handles

%===============================================================================
% Start matrix layout subplotting
%
for j = 1:nChannels
   s = nodesett(j);
   for i = 1:nChannels
      r = nodesett(i);
      if j ~= i || ( j == i && flgPrinting(7) ~= 0)
         if ( j == i && flgPrinting(7) ~= 0) || flgPrinting(4) ~= 0
            h = subplot2(nChannels,nChannels,(i-1)*nChannels+j);
         end
      end
      
%===============================================================================
%     Power Spectra Plots on Main Diagonal
%
      if (j == i) %&& flgPrinting(7) ~= 3 %Power spectrum
         if flgPrinting(7) ~= 0
            SStmp = abs(getCij(SS,r,s,nPlotPoints));
            Ltmp  = abs(getCij(L, r,s,nPlotPoints)); % PDC2/DTF2 estimates (?)
            
            switch flgPrinting(7) % Main diagonal plotting SS and/or PDC2/DTF2
               
               case 1 %Standardized spectra
                  SStmp = SStmp/max(SStmp);
                  h12   = plot(h,w,SStmp);
                  ax(1) = gca;
                  ax(2) = ax(1);
                  if j == 2
                     hh = ylabel('Spectra [a.u.]');
                     pos = get(hh,'Position');
                     if ~isOctave()
                        pos(1) = 0.0055;
                     else
                        pos(1) = (pos(1)+0.00275)/2;
                     end
                     set(hh,'Position',pos, 'FontSize', vFontSize);
                  end
                  set(h12,'LineWidth',vLineWidth,'Color',[0 0 0.7]);

                  
                  if flgScale == 1
                     yTick = [0 .5 1];
                     if i == 1 || i == nChannels
                        yTickLabel = [' 0';'.5';' 1'];
                     else
                        yTickLabel = [' ';' ';' '];
                     end
                     set(h,'XLim', [w_min w_max], 'YLim',kylim, ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                        'YTick',yTick, ...
                        'YTickLabel',yTickLabel, ...
                        'FontWeight','bold',...
                        'FontName', 'Arial', ...
                        'FontSize', vFontSize);
                     if ~isOctave()
                        h.YRuler.TickLabelGapOffset = vTickLabelGapOffset;
                     end
                  else
                     vaxis = double(axis); % axis may return 'single', but oddly
                                           % enough axis() requires 'double'.
                     vaxis =[w_min w_max -0.05*vaxis(4) 1.075*vaxis(4)];
                     set(h,'XLim', [w_min w_max], ...
                        'YLim',[vaxis(3) vaxis(4)], ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                        'FontWeight','bold',...
                        'FontName', 'Arial', ...
                        'FontSize', vFontSize);
                  end
                  
               case 2 % Log spectra on main diagonal
                  SStmp = log(SStmp);
                  SStmp = (SStmp-min(SStmp))/(max(SStmp)-min(SStmp));
                  h12 = plot(h,w,SStmp);
                  ax(1) = gca;
                  ax(2) = ax(1);
    
                  if j == 2
                     hh = ylabel('log Spectra [a.u.]');
                     pos = get(hh,'Position');
                     if ~isOctave()
                        pos(1) = 0.0055;
                     else
                        pos(1) = (pos(1)+0.00275)/2;
                     end
                     set(hh,'Position',pos, ...
                            'FontSize', vFontSize, 'FontWeight', 'normal');
                  end
                  set(h12,'LineWidth',vLineWidth,'Color',[0 0 0.7]);
                  
                  if flgScale == 1
                     yTick = [0 .5 1];
                     if i == 1 || i == nChannels
                        yTickLabel = [' 0';'.5';' 1'];
                     else
                        yTickLabel = [' ';' ';' '];
                     end
                     set(h,'XLim', [w_min w_max], 'YLim',kylim, ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                        'YTick',yTick, ...
                        'YTickLabel',yTickLabel, ...
                        'FontWeight','bold',...
                        'FontName', 'Arial', ...
                        'FontSize', vFontSize);
                     if ~isOctave()
                        h.YRuler.TickLabelGapOffset = vTickLabelGapOffset;
                     end
                  else
                     vaxis = double(axis);
                     vaxis =[w_min w_max -0.05*vaxis(4) 1.075*vaxis(4)];
                     set(h,'XLim', [w_min w_max], ...
                        'YLim',[vaxis(3) vaxis(4)], ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                        'FontWeight','bold',...
                        'FontName', 'Arial', ...
                        'FontSize', vFontSize);
                  end
            end
            
            if j == nChannels
               hxlabel(j) = labelitx(j,chLabels);
            end
            
            set(h,'XTickLabel',[' ';' ';' ';' ';' ';' ']);
            
            if j == 1 && (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
               hylabel(i) = labelity(i,chLabels);
               if nChannels > 10  % Rotate ylabels by 20 degrees if nChannels > 10.
                  set(hylabel(i),'Rotation',70);
               end
            end
         
            ylim = [-0.05 1.075];
            yTick = [0 .5 1];
            if i == nChannels && j == 1
               yTickLabel = [' 0';'.5';' 1'];
            else
               yTickLabel = [' ';' ';' '];
            end                   
            
         else % if flgPrinting(7) == 0 -- main diagonal w/o subplots
            % Put label on last column variable or channel
            if j == nChannels
               hxlabel(j) = labelitx(j,chLabels);
            end
         end % if flgPrinting(7) ~= 0
         
% ==============================================================================      
% p-values plotting - Off diagonal subplots
%
      else % (i~=j) p-values plot

         pvalues_tmp = getCij(pvalues,r,s,nPlotPoints);
         index = pvalues_tmp >= alpha;
         atrib = [0 .4 0]; % Non significant p-values in dark-green line
         kylim  = [-0.05  1.0]; 
         
         hold on
         if flgPrinting(4) == 1 % Linear scale
            plot(h,w,pvalues_tmp, 'Color',atrib,'LineWidth',vLineWidth);
            axis([0 max(w) -0.05  1.0]) 
            hold on
            pvalues_tmp(index) = NaN;
            
            if flgScale
               ylim = kylim;
               yTick = [0 .5 1];
               yTickLabel = [' ';' ';' '];
            end
            
            % Plot significant p-values in red line segments
            if (sum(isnan(pvalues_tmp)) ~= nPlotPoints)
               plot(h,w,pvalues_tmp,'r-','LineWidth',vLineWidth);
            end
            
            if flgPrinting(2) == 1  % Plot significance leve alpha
               if nChannels > 20
                  plot(h,[0 max(w)],[alpha alpha],'k:', ...
                     'LineWidth',vPatnaikLineWidth)
               else
                  plot(h,[0 max(w)],[alpha alpha],'k--', ...
                     'LineWidth',vPatnaikLineWidth)
               end
            end            

         elseif flgPrinting(4) > 1 % log10 scale
            alpha_log10 = log10(alpha);
            pvalues_log10 = log10(abs(pvalues_tmp) + eps);
            plot(h,w,pvalues_log10, 'Color',atrib, 'LineWidth',vLineWidth);
            vaxis = axis;
            if flgScale == 1
               ylim = [-15*1.05 0];
               yTick = [-15 -10 -5 0];
               yTickLabel = [' '; ' '; ' '; ' '];

            elseif flgScale == 2
               if max(abs(pvalues_log10)) <= 5
                  ylim = [-5*1.05 0];
                  yTick = [-5 0];
                  yTickLabel = ['-5'; ' 0'];

               elseif max(abs(pvalues_log10)) <= 10
                  ylim = [-10*1.075 0];
                  yTick = [-10 -5 0];
                  yTickLabel = ['-10'; ' -5'; '  0'];

               else
                  axis([0 max(w) vaxis(3) 0])
                  ylim = [-15*1.075 0];
                  yTick = [-15 -10 -5 0];
                  yTickLabel = ['-15'; '-10'; ' -5'; '  0'];
               end
               
            else
               if max(abs(pvalues_log10)) <= 10
                  vaxis3 = -10;
                  vaxis(3) = 1.075*vaxis3;
               else
                  vaxis3 = -ceil(max(abs(pvalues_log10)));
                  vaxis(3) = 1.075*vaxis3;
               end
               axis([0 max(w) vaxis(3) 0])
               ylim = [vaxis(3) 0];

            end
            hold on
            pvalues_log10(index) = NaN;
            set(h,'YLim',ylim,'YTick', yTick,'YTickLabel',yTickLabel)
            
            if (sum(isnan(pvalues_tmp)) ~= nPlotPoints)
               plot(h,w,pvalues_log10,'r-','LineWidth',vLineWidth);
            end
            if flgPrinting(2) == 1
               if nChannels > 20
                  plot(h,[0 max(w)],[alpha_log10 alpha_log10],'k:', ...
                     'LineWidth',vPatnaikLineWidth)
               else
                  plot(h,[0 max(w)],[alpha_log10 alpha_log10],'k--', ...
                     'LineWidth',vPatnaikLineWidth)
               end
            end
         else
            error('p-values ploting scale value should be either 1 or >1.')
         end
         

         
      end
      
% ==============================================================================
% Labeling axes, and set x- y-axis limits, x- yTickLabels, FontSize 
%
      if (j ~= i && flgPrinting(4) ~= 0) || ( j == i && flgPrinting(7) ~= 0)
         if i == nChannels    % Bottom row subplots
            if j == 1 % subplot with Freq Scale, [bottom left corner subplot]
               if flgScale == 1
                  if flgPrinting(4) == 1
                     set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                        'XTick',2*w_max*[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',{' 0';'  ';'  ';'  ';'  ';str_w_max}, ...
                        'YTick',yTick,'YTickLabel',[' 0'; '.5'; ' 1'], ...
                        'FontSize', vFontSize,'FontWeight','bold')
                     
                  elseif flgPrinting(4) == 2
                     set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                        'XTick',2*w_max*[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',{' 0';'  ';'  ';'  ';'  ';str_w_max}, ...
                        'YTick',yTick,'YTickLabel',['-15'; '-10'; ' -5'; '  0'], ...
                        'FontSize', vFontSize,'FontWeight','bold')
                  end
               else
                  set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                     'XTick',2*w_max*[0 .1 .2 .3 .4 .5], ...
                     'XTickLabel',{' 0';'  ';'  ';'  ';'  ';str_w_max}, ...
                     'FontSize', vFontSize,'FontWeight','bold')                  
               end

               if flgPrinting(4) ~= 0
                  hylabel(i) = labelity(i,chLabels);
               end
               hxlabel(j) = labelitx(j,chLabels);
               
            else % j ~= 1
               if flgScale == 1 %&& flgPrinting(4) == 2 
                  if j ~= nChannels
                     set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                        'YTick',yTick,'YTickLabel',yTickLabel, ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                        'YTick',yTick,'YTickLabel',yTickLabel, ...
                        'FontSize', vFontSize,'FontWeight','bold')
                  end
               else
                  set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                     'XTick',[0 .1 .2 .3 .4 .5], ...
                     'XTickLabel',[' ';' ';' ';' ';' ';' ']);
               end
               set(h,'FontSize', vFontSize,'FontWeight','bold');
               if ~isOctave()
                  h.YRuler.TickLabelGapOffset = vTickLabelGapOffset;
               end
               
               if flgScale == 2
                  set(h,'YTick',yTick,'YTickLabel',yTickLabel, ...
                     'FontSize', vFontSize,'FontWeight','bold', ...
                      'FontName','Arial');
               end
               hxlabel(j) = labelitx(j,chLabels);
            end % j == 1
         elseif i == 1 && j == 2 && flgPrinting(7) == 0 % Special case row #1
            set(h,'XLim', [w_min w_max],'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'FontSize', vFontSize,'FontWeight','bold');
            if flgScale ==2 || flgScale == 1
               set(h,'YTick',yTick,'YTickLabel',yTickLabel, ...
                  'FontSize', vFontSize,'FontWeight','bold', ...
                   'FontName','Arial');
            end
            if (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
               hylabel(i) = labelity(i,chLabels);
               if nChannels > 10  % Rotate ylabels by 20 degrees if numChannels > 10.
                  set(hylabel(i),'Rotation',70);
               end
            end
         elseif j == 1 % Column #1 of lay-out
            set(h,'XLim', [w_min w_max],'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'FontSize', vFontSize,'FontWeight','bold');
            if flgScale == 2, set(h,'YTick',yTick,'YTickLabel',yTickLabel); end
            if i == nChannels
               set(h,'FontSize', vFontSize,'FontWeight','bold', ...
                      'FontName','Arial');
               hxlabel(j) = labelitx(j,chLabels);
            elseif flgScale == 1
             set(h,'XLim', [w_min w_max],'YLim',ylim, ...
                   'YTick',yTick,'YTickLabel',yTickLabel, ...
                   'XTick',[0 .1 .2 .3 .4 .5], ...
                   'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                   'FontSize', vFontSize,'FontWeight','bold');         
            end
            if (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
               hylabel(i) = labelity(i,chLabels);
            end
         else  % j~=1 or j~=nChannels or not(i==1 and j==2 and ~flgPrinting(7))
            set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                  'XTick',[0 .1 .2 .3 .4 .5], ...
                  'XTickLabel',[' ';' ';' ';' ';' ';' ']);

            if flgScale == 1 || flgScale == 2,
               set(h,'YTick',yTick,'YTickLabel',yTickLabel);
            end
            set(h,'FontSize', vFontSize,'FontWeight','bold', ...
                  'FontName','Arial');
         end % i == nChannels
         % Label axes; set x- y-axis limits, x- yTickLabels, FontSize
         
         
% ==============================================================================
         if j == 1
            hylabel(i) = labelity(i,chLabels);
            if nChannels > 20  % Rotate ylabels by 20 degrees if numChannels>20
               set(hylabel(i),'Rotation',70);
            end
         end
         if ~isOctave()
             h.YRuler.TickLabelGapOffset = vTickLabelGapOffset;
         end
% ==============================================================================
% Print GCT p-values right above each subplot
% ==============================================================================
%
         if r ~= s
            if (flgPrinting(4) > 0) && (flgPrinting(5) > 0)
               vaxis = double(axis);
               if isfield(c,'Tragct')
                  if Tragct(r,s)
                     if (flgPrinting(5) > 1)
                        plot(h,0.96*vaxis(2),0.93*vaxis(4)+0.07*vaxis(3),'.', ...
                             'Color',vMarkerColor,'MarkerSize',vMarkerSize);
                     end
                     if isfield(c,'pvaluesgct') && (flgPrinting(5) == 1 || ...
                                                    flgPrinting(5) == 3)
                        if pvaluesgct(r,s) < 1.0e-15
                           text(vaxis(2),1.1*(vaxis(4)-vaxis(3))+vaxis(3), ...
                              ['<1.0e-15'], ...
                              'FontSize',pValuesFontSize, ...
                              'FontWeight','bold','Color',0.6*vMarkerColor, ...
                              'HorizontalAlignment','right')
                        else
                           text(vaxis(2),1.1*(vaxis(4)-vaxis(3))+vaxis(3), ...
                              sprintf('%0.3g',pvaluesgct(r,s)), ...
                              'FontSize',pValuesFontSize,'FontWeight','bold', ...
                              'Color',0.6*vMarkerColor, ...
                              'HorizontalAlignment','right')
                        end
                     end
                  else % not significant GCT
                     if isfield(c,'pvaluesgct') && (flgPrinting(5) == 1 || ...
                                                    flgPrinting(5) == 3)
                        text(vaxis(2),1.1*vaxis(4)-0.1*vaxis(3), ...  % 0.6*vaxis(2)
                             sprintf('%0.3g',pvaluesgct(r,s)), ...
                             'FontSize',pValuesFontSize,'FontWeight','bold', ...
                             'Color',[0 .4 0],'HorizontalAlignment','right')
                     end
                  end
               end
            end  %strrep(num2str(pvaluesgct), '0.', '.')
         end % GCT p-values
      end
   end
end

supAxes = [.08 .075 .84 .84]; % [.12 .11 .84 .80];

%%
% Y-LABEL
%

%[ax3,h3] = suplabel(['log_{10}\textit{p} (\lambda)'],'y',supAxes);
[ax3,h3] = suplabel(['log_{10}p'],'y',supAxes);

switch computer
   case 'PCWIN64'
      set(h3,'FontSize',12,'FontWeight','bold','FontName','Arial')
   case 'GLNXA64'
      set(h3,'FontSize',12,'FontWeight','bold','FontName','Arial')      
   case 'MACI64'
      set(h3,'FontSize',16,'FontWeight','bold','FontName','Arial')
end

pos = get(ax3,'Position');
pos(1) = pos(1) + 0.030;   % 0.0545
set(ax3,'Position',pos);    % Adjust ylabel position

if fs == 1
   strXsuplabel = 'Frequency';
else
   strXsuplabel = 'Frequency (Hz)';
end   

%%
% X-LABEL 
%

[ax1,h1] = suplabel(strXsuplabel,'x',supAxes);
pos = get(ax1,'Position');

switch computer
   case 'PCWIN64'
      set(h1,'FontWeight','bold', 'FontSize',12,'FontName','Arial');
      pos(2) = pos(2) + 0.035;
   case {'GLNXA64','x86_64-pc-linux-gnu','i686-pc-linux-gnu'}
      set(h1,'FontWeight','bold', 'FontSize',12,'FontName','Arial');
      pos(2) = pos(2) + 0.025; %pos(2) + dpos2;           % 0.0545 % 0.025
   case 'MACI64'
      set(h1,'FontWeight','bold', 'FontSize',16,'FontName','Arial');
      pos(2) = pos(2) + 0.010;
   otherwise
      disp('otherwise')
end

if ~isOctave()   
   set(ax1,'Position',pos);    % Adjust xlabel position
else
   set(ax1,'Position',pos);    % Actually suplabel did not work in Octave   
end

% ==============================================================================
% Fine adjustment of axis channel labels positions.
%

for k = 1:nChannels
   set(hxlabel(k),'Units','normalized');
   
   if (flgPrinting(4) ~= 0 || flgPrinting(7) ~= 0)
      set(hylabel(k),'Units','normalized');
      pos = get(hylabel(k),'Position');
      pos(1) = -0.135*nChannels/4;
      set(hylabel(k),'Position',pos);
   end
   
   pos = get(hxlabel(k),'Position');
   pos(2) = -0.145*nChannels/4;
   set(hxlabel(k),'Position',pos);
end

drawnow

end


%% LABELITX 
%          x-axis labeling function
function [hxlabel] = labelitx(j,chLabels) % Labels x-axis plottings
if isOctave()
   vFontSize = 12;
else
   switch computer
      case 'PCWIN64'
         vFontSize = 10;
      case 'GLNXA64'
         vFontSize = 10;
      case 'MACI64'
         vFontSize = 12;
   end
end

if isempty(chLabels)
   hxlabel = xlabel(['j=' int2str(j)]);
   set(hxlabel,'FontSize',vFontSize, 'FontWeight','bold');
   %      'FontName','Arial') % 'FontName','Arial'
else
   hxlabel = xlabel([chLabels{j}]);
   set(hxlabel,'FontSize',vFontSize,'FontWeight','bold');
end

% if isempty(chLabels)
%    if is_octave()
%       hxlabel = xlabel(['j = ' int2str(j)]);
%    else
%       hxlabel = xlabel(['\bf{\it{j}} \rm{\bf{ = ' int2str(j) '}}']);
%    end
%    set(hxlabel,'FontSize',vFontSize,'FontWeight','bold', ...
%         'FontName','Times')
% else
%     hxlabel = xlabel([chLabels{j}]);
%     set(hxlabel,'FontSize',vFontSize,'FontWeight','bold','FontName', 'Arial')
% end

pos = get(hxlabel,'Position');
pos(2) = pos(2) + 0.055;   % 0.0545
set(hxlabel,'Position',pos);
end

%% LABELITY
%   y-axis labeling function
%
function [hylabel] = labelity(i,chLabels) % Labels y-axis plottings
if isOctave()
   vFontSize = 12;
else
   switch computer
      case 'PCWIN64'
         vFontSize = 10;
      case 'GLNXA64'
         vFontSize = 10;
      case 'MACI64'
         vFontSize = 12;
   end
end

if isempty(chLabels)
   hylabel = ylabel(['i=' int2str(i)],...
      'Rotation',90);
   set(hylabel,'FontSize',vFontSize,'FontWeight','bold');
   %      'FontName','Arial')  % 'FontName','Arial', 'Times'
else
   hylabel = ylabel([chLabels{i}]);
   set(hylabel,'FontSize',vFontSize,'FontWeight','bold');
end

pos = get(hylabel,'Position');
pos(1) = pos(1) + 0.0545;   % 0.0545
set(hylabel,'Position',pos)
end
%==========================================================================

