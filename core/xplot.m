%% XPLOT
%        Pretty plot PDC2 / DTF2 (squared-PDC / DTF) connectivity measure
%        and  Coh2 in the frequency domain using matrix layout, optionally with
%        power spectra (and/or auto-PDC2/DTF2) plots along the main diagonal.
%
%% Syntax:
%        [hfigure,hxlabel,hylabel] = XPLOT(strBarTitle,c,flgPrinting,fs,w_max, ...
%                              chLabels,flgColor,flgScale,flgMax,flgSignifColor)
%
%% Input Arguments:
%
%   strBarTitle:  Title for the figure window bar
%
%   c: struct variable with the following fields:
%   |
%   +---.{{pdc|dtf},{pdc2|dtf2},pvalues,th,ci1,ci2,metric,alpha,p,patdenr,
%         patdfr,SS,coh2,pvaluesgct,Tragct}
%
%                 1 2 3 4 5 6 7
%   flgPrinting: [1 1 1 1 1 0 1]; % Example: Plot everything, except coh2.
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
%
%   fs:         sampling frequency
%
%   w_max:      frequency scale upper-limit (<= fs/2, i.e. Nyquist frequency)
%
%   chLabels:   channels id labels (if ==[], channels will be enumerated)
%
%   flgColor:   0: white background;
%               1: white [.1 1], light-blue [.01 .1], purple [0 .01]
%               2: white [.1 1], light-gray [.01 .1], gray [0 .01]
%
%   flgScale:   1: [0 1] / {if max(PDC2/DTF2) > 1}:[0 max(PDC2/DTF2)]
%               2: [0 {if max(PDC/DTF2) > 1}:max(PDC/DTF2)]/[0 1]/[0 .1]/[0 .01]
%                                         based on flgMax (PDC2/DTF2/Thr/CI/all)
%               3: [0 max(PDC2/DTF2/Thr/CI/all)]
%                                         based on flgMax (PDC2/DTF2/Thr/CI/all)
%               4: [0 {max(PDC2/DTF2/Thr/CI/all) or round to {0.01 0.1 1.0}]
%
%   flgMax:     {'PDC'|'DTF'|'Thr'|'CI'|'TCI'|'all'} measure used as upper limit
%               for y-axis scales:
%                    PDC|DTF - PDC2/DTF2 maximum value;
%                        Thr - Patnaik threshold maximum value;
%                        CI  - maximum upper confidence interval value;
%                        TCI - threshold or CI max value;
%                        all - plotted either max([PDC2/DTF2, Thr, CI]) value;
%               See also flgScale.
%
%   flgSignifColor: line color for PDC2/DTF2 plots
%               0: black line
%               1: black / gray  => significant / not PDC2/DTF2
%               2: red   / gray  =>      "
%               3: red   / green =>      "
%               4: red   / black =>      "
%               5: black / green =>      "
%
%% Output:
%     Pretty matrix-layout plots and optionally p-values of GCT
%     hfigure: figure handle
%     hxlabel,hylabel: label's handles
%
%% Description:
%     Cosmetic script to visualize results from 'asymp_pdc.m' or 'asymp_dtf.m'
%     routines in a "standard" matrix-layout.
%
%% Examples:
%   %Try the following example:
%
%   % -------- Example script start here----------------------------------------
%   N = 7;
%   u = randn(2000,N)'; % generate 10 Gaussian random processes.
%   maxIP = 30; alg = 1; criterion = 1; % using Nutall-Strand alg + AIC
%   [IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);
%   nFreqs = 128; metric = 'info'; alpha = 0.05;
%   c = asymp_pdc(u,A,pf,nFreqs,metric,alpha);
%   xplot([],c); shg; % Defaults flgPrinting, fs, w_hmax and flgColor
%   % -------- End here---------------------------------------------------------
%
%   % You might see on average  4 or 5 suplots with red lines, i.e.
%   % significat iPDC at alpha = 0.05, which is the number of analyzed pairs
%   % times alpha, or N(N-1)*alpha = 10(10-1)*0.05 = 4.5 cases or connections
%
%   % Then try the following option:
%   xplot([],c,[1 1 1 0 0 1 1],1,0.5,[],0); shg;
%   % which will plot PDC2/DTF2, threshold, power spectra, coherence plots
%   % with normalized fs=1 Hz; default channel label;
%   % flgColor=0 => no color or rescalling is used.
%
%% See also: XPLOT_TITLE, XPLOT_PVALUES ASYMP_PDC, ASYMP_DTF, SUBPLOT2
  
% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.


function [hfigure,hxlabel,hylabel] = xplot(strBarTitle,c,flgPrinting,fs,w_max, ...
                               chLabels,flgColor,flgScale,flgMax,flgSignifColor)

msgID = 'all';
warning('OFF',msgID);

while 0
   if isOctave()
      disp('Output of ''get(0,''MonitorPositions'') command:')
      get(0, 'MonitorPositions')
      disp('Output of ''get(0,''ScreenSize''')
      get(0, 'ScreenSize')
   else
      disp('1. Output of ''get(0,''MonitorPositions'') command:')
      sz = get(0, 'MonitorPositions');
      [m n] = size(sz);
      if m >= 2,
         fprintf([' a. Multiple-monitor detected.\n'])
         fprintf([' b. Screens resolutions are: \n'])
         disp(sz)
      elseif m == 1
         fprintf([' a. Single monitor detected.\n'])
         fprintf([' b. Screen resolution is: \n'])
         disp(sz)
      else
         error('*** No monitor detected.')
      end

      fprintf('\n2. Output of ''get(0,''ScreenSize''): ')
      sz = get(0, 'ScreenSize');
      disp(sz)
      if m >= 2
         fprintf('\n *** This command can only detect the primary monitor.\n\n')
      end
   end
   pause(5)
   break
end


%%
% <https://www.mathworks.com/matlabcentral/answers/312738-how-to-get-real-screen-size>
% by 
%   <https://www.mathworks.com/matlabcentral/profile/authors/770850>
while 0
   ScreenPixelsPerInch = java.awt.Toolkit.getDefaultToolkit().getScreenResolution()
   ScreenDevices = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().getScreenDevices();
   MainScreen = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getScreen()+1;
   MainBounds = ScreenDevices(MainScreen).getDefaultConfiguration().getBounds();
   MonitorPositions = zeros(numel(ScreenDevices),4);
   for n = 1:numel(ScreenDevices)
      Bounds = ScreenDevices(n).getDefaultConfiguration().getBounds();
      MonitorPositions(n,:) = [Bounds.getLocation().getX() + 1, ...
                    -Bounds.getLocation().getY() + 1 - Bounds.getHeight() + ...
                   MainBounds.getHeight(),Bounds.getWidth(),Bounds.getHeight()];
   end
   disp('MonitorPositions parameter:')
   disp(MonitorPositions)
   pause(5)
   break
end

% ------------------------------------------------------------------------------
% To deal with errobar CapSize option compatibility along different
% MATLAB version R2016b or later (9.1)
if ~isOctave()
   v = ver('MATLAB');
   versionNumberMatlab = str2double(v.Version);
else
   versionNumberMatlab = 8.4; % Does Octave have same behavior as Matlab 2014b?!
end
flgPrint = 'Print';
knargin = 10;     % Number of input arguments
kylim  = [-0.05  1.075]; % y-axis scaling factors at lower and upper limits.
%             kylim1 = [-0.05   1.10]; % Previous values
%             kylim2 = [-0.02   1.05]; % Previous previous values

if nargin < 7, flgColor = 0; end % 0: white background;
if nargin < 8 || isempty(flgScale), flgScale = 1; end
%                        % 1: [0 1] / {if max(PDC2/DTF2) > 1}:[0 max(PDC2/DTF2)]
if nargin < 9, flgMax = 'TCI'; end   % {'PDC'; 'DTF'; 'Thr'; 'CI'; 'TCI'; 'all'}
if nargin < 10, flgSignifColor = 3; end % 3: red/green: signif/~signif PDC2/DTF2
if flgSignifColor < 0 || flgSignifColor > 5, flgSignifColor = 3; end

if nargin <  3,  flgPrinting = [1 0 1 0 0 0 1]; end
if nargin < 4,  fs = 1; end
if nargin <  5,  w_max = fs/2; end

if isfield(c,'dtf2')
   L = c.dtf2;       % |DTF|^2 (nChannels x nChannels x freq)
   if isfield(c,'dtf2_th')
      LTra = c.dtf2_th; % Significant |DTF|^2 on freq range, otherwise NaN.
   end
   flgType = 'DTF';
elseif isfield(c,'pdc2')
   L = c.pdc2;       % |PDC|^2 (nChannels x nChannels x freq)
   if isfield(c,'pdc2_th')
      LTra = c.pdc2_th; % Significant |PDC|^2 on freq range, otherwise NaN.
   end
   flgType = 'PDC';
else
   error('Struct variable c does not hold PDC2/DTF2 analysis results.')
end

if isfield(c,'Tragct')
   Tragct = c.Tragct;
end
if isfield(c,'pvaluesgct')
   pvaluesgct = c.pvaluesgct;
end

SS = c.SS;       % Spectra
Coh2 = c.coh2;   % |Coh|^2
colorCoh2 = [.7 .7 .7];
Lpatnaik = c.th; % Patnaik threshold values for alpha (level of significance)

if isempty(c.ci1) || isempty(c.ci2)
   flgPrinting(4) = 0;
else
   L2vlower = c.ci1;
   L2vupper = c.ci2;
end

metric = c.metric;      % Valid options: "euc", "diag" or "info"
[nChannels,~,nFreqs] = size(L); % nChannels = ~ dummy = # of channels/time series;
%                       % nFreqs = # of points calculated on frequency scale

nodesett = 1:nChannels;

% ------------------------------------------------------------------------------
% figure scaling parameter adjustment
vChannels = max([10, nChannels]);

% Formatting parameters for p-values plotting above each subplot.
N1 = 1; N30 = 30; V1 = 25; V30 = 10; % was 10,5
vMarkerSize = V1 + (vChannels - N1)*(V30 - V1)/(N30 - N1);

vMarkerColor = [ 1       0         1]; % Purple
N1 = 10; N30 = 24; V1 = 8; V30 = 5; % was 10,5
pValuesFontSize = V1 + (vChannels - N1)*(V30 - V1)/(N30 - N1);

% if nChannels >20
%    pValuesFontSize = 4;
% elseif nChannels > 15
%    pValuesFontSize = 5;
% else
%    pValuesFontSize = 7;
% end

if w_max > fs/2 + eps
   error(['The parameter w_max should be <= Nyquist frequency,' ...
                                                        'i.e, w_max <= fs/2.'])
end

% ==============================================================================
%          Check flgPrinting choices and availabitiy of proper measures
% ==============================================================================
%   flgPrinting: [X X X X X X X];
%       blue-line | | | | | | 7--Spectra (0: w/o; 1: Linear; 2: Log; 3: PDC2;
%                 | | | | | |      4: Linear normalized; 5: Log spectra + PDC2)
%            gray | | | | | 6--Coh2 (0: w/o Coh2; 1: with Coh2)
%     dark-purple | | | | 5--Mark GCT significant pairs and p-values printing
%                 | | | |      (0:w/o; 1: p-values print; 2: dot plot +GCT;
%                 | | | |      3: p-values print + dot plot for significant GCT)
%     dashed-blue | | | 4--Confidence interval (0:w/o; 1: Dashed-lines;
%                 | | |          2: Shaded-plot; 3: Error-bar 1; 4: Error-bar 2
%            red  | | 3--{0:1} Significant PDC2|DTF2 in red lines
%    dashed-black | 2--{0:1} Patnaik threshold level in black dashed-lines
%           green 1-- {0:1} PDC2/DTF2 in green lines or black w/o statistics.

if (flgPrinting(1) == 1) && (flgPrinting(5) == 1) && ~isfield(c,'Tragct')
   fprintf(2,['\n*** Warning: to print GCT p-values, assign ''gct_alg()''' ...
                                               'results to c struct. ***\n'])
   fprintf(2,['\n    Tr_gct, pValue_gct, Tr_igct, pValue_igct] =' ...
                               'gct_alg(u,A,pf,gct_signif,flgPrintResults);'])
   fprintf(2,['\n   then, before calling xplot(), execute:'])
   fprintf(2,['\n    c.Tragct = Tr_gct;'])
   fprintf(2,['\n    c.pvaluesgct = pValue_gct;\n\n'])
end

if sum(flgPrinting([1 6 7])) == 0
   error(['''xplot.m'' requires at least one measure to plot: PDC/DTF, ' ...
          'Coh2 and/or SS. flgPrinting=[' int2str(flgPrinting) ']'])
end

if c.alpha <= 0
   if isfield(c,'pdc2') || isfield(c,'dtf2')
      if flgPrinting(1)
         flgPrinting = [1 0 0 0 0 1 1] .* flgPrinting;
      else
         flgPrinting(6) = 1;
         flgPrinting = [0 0 0 0 0 1 1] .* flgPrinting;
      end
      L2vupper = [];
      L2vlower = [];
   else
      error('No PDC2 or DTF2 measure available in result struct variable.')
   end
else  % alpha > 0. For fast routine asymp_pdc routine and flgPrinting(4)==0
   if isfield(c,'pdc2')
      if isempty(c.ci1)
         flgPrinting = [1 1 1 0 1 1 1] .* flgPrinting;
         L2vupper = [];
         L2vlower = [];
      end
   elseif isfield(c,'dtf2')
      %nop
   else
      error('No PDC2 or DTF2 measures available in result struct variable.')
   end
   if flgPrinting(1) == 0
      flgPrinting = flgPrinting .* [0 0 0 0 0 1 1];
      if flgPrinting(7) == 4, flgPrinting(7) = 1; end
      if flgPrinting(7) == 5, flgPrinting(7) = 2; end
      if flgPrinting(7) >  5, flgPrinting(7) = 1; end
      if flgPrinting(7) == 3, flgPrinting(7) = 1; end
   end
end

if exist('shadedplot.m','file') ~= 2, flgColor = 0; end


%%
% Figure setting
%

hfigure = figure;

set(hfigure,'PaperOrientation','landscape', ...
            'renderer','painters'); % To get vector figure

if isempty(strBarTitle)
   set(hfigure,'NumberTitle','off','MenuBar','none');
else
   set(hfigure,'NumberTitle','off','MenuBar','none','Name', strBarTitle);
end

w = 0:fs/(2*nFreqs):w_max-fs/(2*nFreqs);
nPlotPoints = length(w);
w_min = w(1);

% ------------------------------------------------------------------------------
% % adjust ticklabels away from axes
% https://www.mathworks.com/matlabcentral/answers/2318-set-position-of-tick-labels
%
% a=gca;
% a.XRuler.TickLabelGapOffset = -8;
% % negative numbers move the ticklabels down (positive -> up)
% a.YRuler.TickLabelGapOffset = -8;
% % negative numbers move the ticklabels right (negative -> left)

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
         vTickLabelGapOffset = -1.5; % YTickLabel distance from y-axis  =1.5
      case {'GLNXA64','x86_64-pc-linux-gnu','i686-pc-linux-gnu'}
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

% ==============================================================================
%               FontSize scaling parameter adjustment
% ==============================================================================

vFontSize = (vChannels - N0)*(f1 - f0)/(N1 - N0) + f0; % YTickLabel Font Size

if vChannels >15
   vLineWidth = 1.25;
   vPatnaikLineWidth = 0.75;
else
   vLineWidth = 1.75;
   vPatnaikLineWidth = 1.25;
end

% ==============================================================================
%              Set figures size for subplot layout
% ==============================================================================

% Screen dimension in pixel.
% Read the following:
% <https://undocumentedmatlab.com/articles/working-with-non-standard-dpi-displays>

set(0,'units','pixels');
sz = get(0,'ScreenSize');

skreenfactor = 2.4; % Ad hoc factor that should 

if ~isOctave()
%    https://www.mathworks.com/matlabcentral/answers/312738-how-to-get-real-screen-size
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

%%
% Let's adopt printing paper size proportion as A4
A4width = 29.7; A4height = 21.0; % A4 landscape paper orientation size in centimeterw
A4ratio = sqrt(2);

% Real figure size for printing or publication in centimeters
width = min([A4width * vChannels/25  A4width]);
height = min([A4height * vChannels/25 A4height]);

% What follow is a kludge solution to determine figure size in normalized units
% that might allow handling the cases of multiple monitors set up in Octave.
% % % rwidth  = pxwidthscreen/sz(3)/khmon/kvmon;
% % % rheight = pxheightscreen/sz(4)/khmon/kvmon*screenratio1/screenratio0;

% rwidth = rwidth * width/A4width;
% rheight = rheight * height/A4height - 0.03;

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
         'position',[(1/khmon-rwidth)/2 (1-1.087*rheight) rwidth rheight]); % -0.03

   case {'GLNXA64','x86_64-pc-linux-gnu','i686-pc-linux-gnu'}
      set(hfigure,'units','normalized', ...
         'position',[(1/khmon-rwidth)/2 (1-1.087*rheight) rwidth rheight]);

   case 'MACI64'
      set(hfigure,'units','normalized', ...
         'position',[(1/khmon-rwidth)/2 (1-1.087*rheight) rwidth rheight]);

   otherwise
      error('computer function not working properly.')
end

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
   str_w_max = sprintf('%5.0f', w_max);
else
   str_w_max = 'fs/2';
end

if nargin < (knargin-4) || isempty(chLabels)
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
      chLabels = c.chLabels;
   else
      chLabels = [];
   end
   if ~isempty(chLabels) && max(size(chLabels)) < nChannels
      error('2 NOT ENOUGH CHANNEL LABELS 2.');
   else
      disp('3 NOT ENOUGH CHANNEL LABELS. Default labels assumed.');
      % Default labeling consist in enumerating the channels from 1 to nChannels
   end
end

hxlabel = 0; % x-axis labels' handles
hylabel = 0; % y-axis labels' handles


if c.alpha == 0     % With no statistics, PDC2 is
   atrib1 = [0 0 0]; % plotted in black lines.
   atrib2 = [0 0 0];
else
   grayline = 0.71;
   switch flgSignifColor
      case 0
         atrib1 = [0 0 0];  %black line
         atrib2 = [0 0 0];  %black line
      case 1
         atrib1 = [0 0 0];  % 1: black / gray -> signif /not signif PDC2/DTF2
         atrib2 = grayline*[1 1 1];  %gray line
      case 2
         atrib1 = [1 0 0];  % 2: red  /  gray ->      "
         atrib2 = grayline*[1 1 1];  %gray line
      case 3
         atrib1 = [1 0 0];  % 3: red  / green ->      "
         atrib2 = [0 1 0];  %    green ->      "
      case 4
         atrib1 = [1 0 0];  % 4: red        "
         atrib2 = [0 0 0];  %    black ->      "
      otherwise
         atrib1 = [0 0 0];  % 5: black        "
         atrib2 = [0 .5 0];  %    green ->      "
   end
end

% ==============================================================================
%     Loop for nChannels x nChannels Lay-out DTF2 or PDC2 plots in the
%                        frequency domain
% ==============================================================================
%
for j = 1:nChannels
   s = nodesett(j);
   for i = 1:nChannels
      r = nodesett(i);
      if j ~= i || ( j == i && flgPrinting(7) ~= 0)
         if ( j == i && flgPrinting(7) ~= 0) || flgPrinting(1) ~= 0 ...
                                             || flgPrinting(6) ~= 0
            h = subplot2(nChannels,nChannels,(i-1)*nChannels+j);
         end
      end

% ==============================================================================
%               Power Spectrum Plotting on Main Diagonal
% ==============================================================================

      if (j == i) %&& flgPrinting(7) ~= 3 %Power spectrum
         if flgPrinting(7) ~= 0
            SStmp = abs(getCij(SS,r,s,nPlotPoints));
            Ltmp  = abs(getCij(L, r,s,nPlotPoints)); % PDC2/DTF2 estimates (?)

            switch flgPrinting(7) % Main diagonal plotting SS and/or PDC2/DTF2

               case 1 %Standardized spectra
                  SStmp = SStmp/max(SStmp);
                  h12   = plot(w,SStmp);
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
                     set(hh,'Position',pos, 'FontSize', vFontSize, ...
                            'FontName','Arial');
                  end
                  set(h12,'LineWidth',vLineWidth,'Color',[0 0 0.7]);


                  if flgScale == 1 || flgScale == 2
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
                        'FontName', 'Arial', ... %'Helvetica-Narrow', ...
                        'FontSize', vFontSize);
                     if ~isOctave()
                        h.YRuler.TickLabelGapOffset = vTickLabelGapOffset;
                     end
                  else
                     vaxis = double(axis); % axis may return single, and oddly
                                           % enough axis() requires doubles.
                     vaxis =[w_min w_max -0.05*vaxis(4) 1.10*vaxis(4)];
                     set(h,'XLim', [w_min w_max], ...
                        'YLim',[vaxis(3) vaxis(4)], ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                        'FontWeight','bold',...
                        'FontName', 'Arial', ... %'Helvetica-Narrow', ...
                        'FontSize', vFontSize);
                  end

               case 2 % Log spectra on diagonal
                  SStmp = log(SStmp);
                  SStmp = (SStmp-min(SStmp))/(max(SStmp)-min(SStmp));
                  h12 = plot(w,SStmp);
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
                            'FontSize', vFontSize, 'FontWeight', 'normal', ...
                            'FontName', 'Arial') %'Helvetica-Narrow', ...);
                  end
                  set(h12,'LineWidth',vLineWidth,'Color',[0 0 0.7]);

                  if flgScale == 1 || flgScale == 2
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
                        'FontName',  'Arial', ... %'Helvetica-Narrow', ...
                        'FontSize', vFontSize);
                     if ~isOctave()
                        h.YRuler.TickLabelGapOffset = vTickLabelGapOffset;
                     end
                  else
                     vaxis = double(axis);
                     vaxis =[w_min w_max -0.05*vaxis(4) 1.10*vaxis(4)];
                     set(h,'XLim', [w_min w_max], ...
                        'YLim',[vaxis(3) vaxis(4)], ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                        'FontWeight','bold',...
                        'FontName', 'Arial', ... %'Helvetica-Narrow', ...
                        'FontSize', vFontSize);
                  end

               case 3 % PDC2/DTF2 on main diagonal
                  if flgPrinting(1)
                     Ltmp = abs(getCij(L,r,s,nPlotPoints)); % PDC2/DTF2 values
                     plot(w,Ltmp,'Color',atrib2,'LineWidth',vLineWidth);

                     hold on
                     ax(1) = gca;
                     ax(2) = ax(1);

                     if c.alpha ~= 0
                        % Significant PDC2/DTF2
                        LTratmp = abs(getCij(LTra,r,s,nPlotPoints));
                        if flgPrinting(2) == 1
                           atrib='k--'; % Patnaik threshold in black dashed line
                           plot(w,abs(getCij(Lpatnaik,r,s,nPlotPoints)), ...
                                           atrib,'LineWidth',vPatnaikLineWidth);
                        end

                        if flgPrinting(3)
                           if flgPrinting(4)
                              if isempty(c.ci1)
                                 L2vuppertmp   = Ltmp;
                              else
                                 L2vuppertmp = abs(getCij(L2vupper,r,s, ...
                                                                  nPlotPoints));
                              end
                              flgSignif      = sum(~isnan(LTratmp));
                              indexSignif    = ~isnan(LTratmp);
                              indexNotSignif = isnan(LTratmp);
                              % CI plotting on main diagonal for PDC or DTF
                              plotCI(flgPrinting,flgSignif,w,r,s,Ltmp, ...
                                           L2vupper,L2vlower,nPlotPoints, ...
                                           indexNotSignif,vPatnaikLineWidth, ...
                                           versionNumberMatlab);
                           end % flgPrinting(4) CI

                           % Draw significant PDC2/DTF2 in different colors
                           plot(h,w,LTratmp,'Color',atrib1, ...
                                                        'LineWidth',vLineWidth);

                           if flgScale == 1 || flgScale == 2 || flgScale == 3
                              yTick = [0 .5 1];
                              if i == 1 || i == nChannels
                                 yTickLabel = [' 0';'.5';' 1'];
                              else
                                 yTickLabel = [' ';' ';' '];
                              end
                              set(h,'XLim',[w_min w_max],'YLim',kylim, ...
                                 'XTick',[0 .1 .2 .3 .4 .5], ...
                                 'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                                 'YTick',yTick, ...
                                 'YTickLabel',yTickLabel, ...
                                 'FontWeight','bold',...
                                 'FontName', 'Arial', ... %'Helvetica-Narrow', ...
                                 'FontSize', vFontSize);
                              if ~isOctave()
                                 h.YRuler.TickLabelGapOffset = ...
                                                            vTickLabelGapOffset;
                              end
                           else
                              vaxis = double(axis);
                              vaxis =[w_min w_max -0.05*vaxis(4) 1.10*vaxis(4)];
                              set(h,'XLim', [w_min w_max], ...
                                 'YLim',[vaxis(3) vaxis(4)], ...
                                 'XTick',[0 .1 .2 .3 .4 .5], ...
                                 'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                                 'FontWeight','bold',...
                                 'FontName', 'Arial', ... %'Helvetica-Narrow', ...
                                 'FontSize', vFontSize);
                           end
                        end % flgPrinting(3) Significant PDC2/DTF2
                     end  % c.alpha ~= 0

                     plot(h,w,Ltmp,'Color',atrib2,'LineWidth',vLineWidth);
                     if c.alpha ~= 0
                        plot(h,w,LTratmp,'Color',atrib1,'LineWidth',vLineWidth);
                     end

                     if j == nChannels
                        hxlabel(j) = labelitx(j,chLabels);
                     end

                     set(h,'XTickLabel',[' ';' ';' ';' ';' ';' ']);

                     if j == 1 && (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
                        hylabel(i) = labelity(i,chLabels);
                     end

                     set(h,'FontSize', vFontSize,'FontWeight','bold', ...
                           'FontName', 'Arial') %, 'Helvetica-Narrow');
                     if ~isOctave()
                        h.YRuler.TickLabelGapOffset = vTickLabelGapOffset;
                     end
                  else

                  end

               case {4,5} % Log spectra + PDC2/DTF2
                  Ltmp = abs(getCij(L,r,s,nPlotPoints));
                  maxPDC = max([1 max(Ltmp)]);

                  % Power spectra
                  if flgPrinting(7) == 4     % Spectra + PDC2/DTF2
                     SStmp = SStmp/max(SStmp);
                  else % flgPrinting(7) == 5 % Log spectra + PDC2/DTF2
                     SStmp = log(SStmp);
                     SStmp = (SStmp-min(SStmp))/(max(SStmp)-min(SStmp));
                  end

                  [ax,h1,h2] = plotyy(h,w,SStmp,w,Ltmp);
                  set(h2,'Color',atrib2,'LineWidth',vLineWidth);

                  hold on

                  if c.alpha ~= 0
                     % Significant PDC2/DTF2
                     LTratmp = abs(getCij(LTra,r,s,nPlotPoints));
                     if flgPrinting(2) == 1
                        atrib='k--'; % Patnaik threshold in black dashed line
                        plot(w,abs(getCij(Lpatnaik,r,s,nPlotPoints)), ...
                           atrib,'LineWidth',vPatnaikLineWidth);
                     end

                     if flgPrinting(3)
                        if flgPrinting(4)
                           if isempty(c.ci1)
                              L2vuppertmp = Ltmp;
                           else
                              L2vuppertmp = abs(getCij(L2vupper,r,s, ...
                                                                 nPlotPoints));
                           end
                           flgSignif      = sum(~isnan(LTratmp));
                           indexSignif    = ~isnan(LTratmp);
                           indexNotSignif = isnan(LTratmp);

                           % CI plotting on main diagonal for PDC or DTF
                           plotCI(flgPrinting,flgSignif,w,r,s,Ltmp, ...
                                           L2vupper,L2vlower,nPlotPoints, ...
                                           indexNotSignif,vPatnaikLineWidth, ...
                                           versionNumberMatlab);

                        end % flgPrinting(4) CI

                        if flgScale == 1
                           yTick = [0 .5 1];
                           if i == 1 || i == nChannels
                              yTickLabel = [' 0';'.5';' 1'];
                           else
                              yTickLabel = [' ';' ';' '];
                           end
                           set(h,'YTick',yTick,'YTickLabel',yTickLabel);
                        end
                     end % flgPrinting(3) Significant PDC2/DTF2
                  end  % c.alpha ~= 0

                  if j == 2
                     if flgPrinting(7) == 4
                        hh = ylabel('Spectra [a.u.]');
                     else
                        hh = ylabel('log Spectra [a.u.]');
                     end
                     pos = get(hh,'Position');
                     if ~isOctave()
                        pos(1) = 0.005;
                     else
                        pos(1) = (pos(1)+0.0025)/2;
                     end
                     set(hh,'Position',pos,'FontSize', vFontSize);
                  end

                  ylim  = maxPDC * kylim;

                  set(h1,'LineWidth',vLineWidth,'Color',[0  0 0.7]);
                  set(h2,'LineWidth',vLineWidth,'Color',[0  1   0]);
                  yTick = [0 .5 1.0];
                  yTickLabel = [' ';' ';' '];
                  if j == nChannels, yTickLabel = [' 0';'.5';' 1']; end
                  set(ax(2),'XLim', [w_min w_max], 'YLim',ylim, ...
                     'XTick',[0 .1 .2 .3 .4 .5], ...
                     'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                     'YTick',yTick, ...
                     'YTickLabel',yTickLabel, ...
                     'FontWeight','bold',...
                     'FontName','Arial', ... %'FontName', 'Helvetica-Narrow', ...
                     'FontSize', vFontSize);
                  if j == 1, yTickLabel = [' 0';'.5';' 1']; end
                  set(ax(1), 'XLim', [w_min w_max], 'YLim',ylim, ...
                     'YTick',yTick, ...
                     'YTickLabel',yTickLabel, ...
                     'FontWeight','bold',...
                     'FontName','Arial', ... % 'Helvetica-Narrow', ...
                     'FontSize', vFontSize);
                  if ~isOctave()
                     ax(1).YRuler.TickLabelGapOffset = vTickLabelGapOffset;
                     ax(2).YRuler.TickLabelGapOffset = vTickLabelGapOffset;
                  end

                  hold(ax(1), 'all'); hold(ax(2), 'all');

                  plot(h,w,Ltmp,'Color',atrib2,'LineWidth',vLineWidth);
                  if c.alpha ~= 0
                     plot(ax(2),w,LTratmp,'Color',atrib1,'LineWidth',vLineWidth);
                  end


                  hold on
                  %            grid on % For plots along the main diagonal.
                  ylim  = kylim;

                  set(gca,'Box'         , 'on'     , ...
                          'TickDir'     , 'in'     , ...
                          'TickLength'  , [.02 .02], ...
                          'XMinorTick'  , 'off'    , ...
                          'YMinorTick'  , 'off');

                  set(ax(1),'XLim',[w_min w_max], ...
                            'XTick',[0 .1 .2 .3 .4 .5],'XTickLabel',[' '],...
                            'YLim',ylim); %,'YTick',[]);
                  if j == nChannels && (flgPrinting(7) == 4 || ...
                                                            flgPrinting(7) == 5)
                     set(ax(2),'XLim',[w_min w_max],'XTick',[0 .1 .2 .3 .4 .5],...
                        'XTickLabel',[' '],'YLim', ylim);
                  end

                  set(h, 'Color', 0.9*[1 1 1],'layer','bottom'); % Background color
                  set(h, 'TickLength', [.02 .02]);
                  if j == nChannels
                     hxlabel(j) = labelitx(j,chLabels);
                  end

                  set(h,'XTickLabel',[' ';' ';' ';' ';' ';' ']);

                  if j == 1 && (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
                     hylabel(i) = labelity(i,chLabels);
                  end

                  if c.alpha ~= 0
                     plot(w,LTratmp,'Color',atrib1,'LineWidth',vLineWidth);
                  end


                  if flgScale == 1 || flgScale == 2  % *****
                     yTick = [0 .5 1];
                     if i == 1 || i == nChannels
                        yTickLabel = [' 0';'.5';' 1'];
                     else
                        yTickLabel = [' ';' ';' '];
                     end
                     set(h,'YTick',yTick,'YTickLabel',yTickLabel);
                  end
            end

            if j == nChannels
               hxlabel(j) = labelitx(j,chLabels);
            end

            set(h,'XTickLabel',[' ';' ';' ';' ';' ';' ']);

            if j == 1 && (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
               hylabel(i) = labelity(i,chLabels);
            end

            ylim = [-0.05 1.075];

            yTick = [0 .5 1];
            if i == nChannels && j == 1
               yTickLabel = [' 0';'.5';' 1'];
            else
               yTickLabel = [' ';' ';' '];
            end

         else % if flgPrinting(7) == 0 Main Diagonal w/o subplots
            % Label last column variable or channel
            if j == nChannels
               hxlabel(j) = labelitx(j,chLabels);
            end
         end % if flgPrinting(7) ~= 0

% ==============================================================================
%     PDC2 and Coh2 plotting - Off diagonal subplots
% ==============================================================================

      else % (i~=j) PDC2 and Coh2  % if NOT {(j == i) %&& flgPrinting(7) ~= 3}
         if flgPrinting(6)
            Coh2tmp = abs(getCij(Coh2,r,s,nPlotPoints));
            % Omnibus options
            ylim = kylim;
            yTick = [0 .5 1];
            yTickLabel = ['0';' ';'1'];
         end
         if flgPrinting(1)
            Ltmp    = abs(getCij(L,r,s,nPlotPoints));
            if c.alpha ~= 0
               LTratmp = abs(getCij(LTra,r,s,nPlotPoints));
               Lpatmp  = abs(getCij(Lpatnaik,r,s,nPlotPoints));
               if isempty(c.ci1)
                  L2vuppertmp   = Ltmp;
               else
                  L2vuppertmp = abs(getCij(L2vupper,r,s,nPlotPoints));
               end
               flgSignif    = sum(~isnan(LTratmp));
               indexSignif  = ~isnan(LTratmp);
               indexNotSignif = isnan(LTratmp);
               atrib2 = [0.41,0.41,0.41];
            else % c.alpha == 0     % With no statistics, PDC2 is
               atrib1 = [0 0 0]; % plotted in black lines.
               atrib2 = [0 0 0];
               flgSignif = 0;
               flgSignifColor = 0;
            end
            grayline = 0.71;

            switch flgSignifColor
               case 0
                  atrib1 = [0 0 0];  %black line
                  atrib2 = [0 0 0];  %black line
               case 1
                  atrib1 = [0 0 0];  % 1: black / gray -> significant
                                     %                     /not signif PDC2/DTF2
                  atrib2 = grayline*[1 1 1];  %gray line
               case 2
                  atrib1 = [1 0 0];  % 2: red  /  gray ->      "
                  atrib2 = grayline*[1 1 1];  %gray line
               case 3
                  atrib1 = [1 0 0];  % 3: red  / green ->      "
                  atrib2 = [0 1 0];  %    green ->      "
               case 4
                  atrib1 = [1 0 0];  % 4: red        "
                  atrib2 = [0 0 0];  %    black ->      "
               otherwise
                  atrib1 = [0 0 0];  % 5: black        "
                  atrib2 = [0 1 0];  %    green ->      "
            end % switch flgSignifColor

            plot(h,w,Ltmp,'Color',atrib2,'LineWidth', vLineWidth);
            ax(1) = gca;
            ax(2) = ax(1);
            grid off
            hold on

            if c.alpha ~= 0
               if isempty(c.ci1)
                  maxCI = max(Ltmp);
               else
                  maxCI  = max(L2vuppertmp(indexSignif));
               end
               maxPDC = max(Ltmp);
               maxThr = max(Lpatmp);
               if flgPrinting(6) ~= 0
                  maxCoh = max(Coh2tmp);
               end
            else
               maxCI  = max(Ltmp);
               maxPDC = maxCI;
               maxThr = maxCI;
               if flgPrinting(6) ~= 0
                  maxCoh = max(Coh2tmp);
               end
            end
% ==============================================================================
%           First flgScale switch
% ==============================================================================
%
            switch flgScale
               case 1 % [0 1] / [0 max(PDC2)]
                  if maxPDC <= 1
                     ylim = kylim;
                     maxValue = 1.0;
                  else
                     ylim = maxPDC*kylim;
                     maxValue = maxPDC;
                  end

               case {0,2,3,4} % [0 max(PDC2)]/[0 1]/[0 .1]/[0 .01]
                  %ytickformat('.%02d');
                  switch upper(flgMax)
                     case 'PDC'
                        maxValue = maxPDC;

                     case 'THR'
                        if flgPrinting(2)
                           maxValue = max([maxThr maxPDC]);
                        else
                           maxValue = maxPDC;
                        end

                     case 'CI'
                        if flgSignif && flgPrinting(4)
                           maxValue = max([maxCI maxPDC]);
                        else
                           maxValue = maxPDC;
                        end

                     case 'COH'
                        if flgPrinting(6)
                           maxValue = max([maxCoh maxPDC]);
                        else
                           maxValue = maxPDC;
                        end

                     case 'TCI'
                        if (flgSignif && flgPrinting(4)) && flgPrinting(2)
                           maxValue = max([maxPDC maxThr maxCI]);
                        elseif flgSignif && flgPrinting(4)
                           maxValue = max([maxPDC maxCI]);
                        elseif flgPrinting(2)
                           maxValue = max([maxPDC maxThr]);
                        else
                           maxValue = maxPDC;
                        end

                     case 'ALL'
                        if (flgSignif && flgPrinting(4)) && flgPrinting(2) ...
                                                               && flgPrinting(6)
                           maxValue = max([maxPDC maxThr maxCI maxCoh]);
                        elseif (flgSignif && flgPrinting(4)) && flgPrinting(2)
                           maxValue = max([maxPDC maxThr maxCI]);
                        elseif (flgSignif && flgPrinting(4)) && flgPrinting(6)
                           maxValue = max([maxPDC maxCoh maxCI]);
                        elseif flgPrinting(2) && flgPrinting(6)
                           maxValue = max([maxPDC maxThr maxCoh]);
                        elseif flgPrinting(2)
                           maxValue = max([maxPDC maxThr]);
                        elseif (flgSignif && flgPrinting(4))
                           maxValue = max([maxPDC maxCI]);
                        elseif flgPrinting(6)
                           maxValue = max([maxPDC maxCoh]);
                        else
                           maxValue = maxPDC;
                        end
                     otherwise
                        maxValue = 1.0;
                  end  % switch upper(flgMax)
            end % switch flgScale


% ==============================================================================
%           Second flgScale switch
% ==============================================================================
%
            switch flgScale
               case 1
                  yTick = [0 .5 1];
                  ylim = kylim;
                  if i == nChannels && j == 1
                     yTickLabel = [' 0';'.5';' 1'];
                  else
                     yTickLabel = [' ';' ';' '];
                  end
               case 2
                  if maxValue < 0.0001
                     ylim = 0.002*kylim; % [-0.0001    0.0022];
                     yTick = [0 .001];
                     yTickLabel = ['   0';'.001'];
                  elseif maxValue < 0.0105
                     ylim = .01*kylim; % [-0.0005    0.0110
                     yTick = [0 .005 .01];
                     yTickLabel = ['  0';'   ';'.01'];
                  elseif maxValue < 0.105
                     ylim = .1*kylim; % [-0.005 0.11];
                     yTick = [0 .05 .1];
                     yTickLabel = [' 0';'  ';'.1'];
                  elseif maxValue < 1.05
                     ylim = kylim;
                     yTick = [0 .5 1];
                     yTickLabel = ['0';' ';'1'];
                  else
                     ylim = maxValue*kylim;
                     yTick = [0 .5 1];
                     yTickLabel = ['0';' ';'1'];
                  end
               case 3
                  ylim = maxValue * kylim;
                     yTick = [0 .5 1];
                     yTickLabel = ['0';' ';'1'];
               case 4
                  if maxValue < 0.0001
                     ylim = 0.001*kylim; %[-0.00004 0.00210];
                  elseif maxValue < 0.005
                     ylim = maxValue*kylim; % [-0.02 1.05];
                  elseif maxValue < 0.0105
                     ylim = 0.01*kylim; %[-0.0005 0.0110]; %[-0.05  1.075]
                  elseif maxValue < 0.05
                     ylim = maxValue*kylim;
                  elseif maxValue < 0.105
                     ylim = 0.1*kylim; %[-0.005 0.110];
                  elseif maxValue < 0.5
                     ylim = maxValue*kylim;
                  elseif maxValue < 1.05
                     ylim = kylim;
                  else
                     ylim = maxValue*kylim;
                  end
            end
% ==============================================================================
% flgColor --  Set subplots background color
% ==============================================================================
%
            set(gca, 'Layer', 'bottom');
            switch flgColor
               case 0
                  % White background.

               case 1
                  if maxValue < 0.001
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        [0.7 0.7 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.8*[1 0.7 1]);
                  elseif maxValue < 0.01
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        [0.7 0.7 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.8*[1 0.7 1]);
                  else
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        [0.7 0.7 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.8*[1 0.7 1]);
                  end

               case 2
                  if maxValue < 0.001
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        0.9*[1 1 1]);  %light gray
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.75*[1 1 1]);  % darker gray
                  elseif maxValue < 0.01
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        0.9*[1 1 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.75*[1 1 1]);
                  else
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        0.9*[1 1 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.75*[1 1 1]);
                  end
            end % switch flgColor



         else % flgPrinting(1) == 0
%
%             disp('In test 2 ...')

            if exist('ylim') ~= 1,
               ylim = kylim;
            end
            yTick = [0 .5 1];
            if i == nChannels && j == 1
               yTickLabel = [' 0';'.5';' 1'];
            else
               yTickLabel = [' ';' ';' '];
            end
         end % if flgPrinting(1)

         grid off

%===============================================================================
% Labeling axes, and set x- y-axis limits, x- yTickLabels, FontSize
% ==============================================================================
%
         if i == nChannels    % Bottom row subplots
            if j == 1 % subplot with Freq Scale
               if flgColor == 4
                  if max(Ltmp) < 0.001  % Review, else is identical
                     set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' 0';'  ';'  ';'  ';'  ';str_w_max], ...
                        'FontSize', vFontSize,'FontWeight','bold', ...
                        'FontName','Arial');
                     if flgScale == 2
                        set(h,'YTick',yTick,'YTickLabel',yTickLabel, ...
                           'FontSize', vFontSize,'FontWeight','bold', ...
                           'FontName','Arial') %'Helvetica-Narrow');
                     end
                  else
                     set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' 0';'  ';'  ';'  ';'  ';str_w_max], ...
                        'FontSize', vFontSize,'FontWeight','bold');
                     if flgScale ==2
                        set(h,'YTick',yTick,'YTickLabel',yTickLabel, ...
                           'FontSize', vFontSize,'FontWeight','bold', ...
                           'FontName', 'Arial') %'Helvetica-Narrow');
                     end
                  end % max(Ltmp) < 0.001 Review
               else
                  set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                     'XTick',2*w_max*[0 .1 .2 .3 .4 .5], ...
                     'XTickLabel',{' 0';'  ';'  ';'  ';'  ';str_w_max}, ...
                     'FontSize', vFontSize,'FontWeight','bold', ...
                     'FontName','Arial');
                  if flgScale ==2 || flgScale == 1
                     set(h,'YTick',yTick,'YTickLabel',yTickLabel, ...
                        'FontSize', vFontSize,'FontWeight','bold', ...
                        'FontName','Arial') %'Helvetica-Narrow');
                  end
               end  % flgColor == 4
               if (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
                  hylabel(i) = labelity(i,chLabels);
               end
               hxlabel(j) = labelitx(j,chLabels);
            else % j ~= 1
               if flgScale == 1
                  if j ~= nChannels
                     set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                        'YTick',yTick,'YTickLabel',yTickLabel, ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' ';' ';' ';' ';' ';' ']);
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
                     'FontName','Arial') %'Helvetica-Narrow');
               end
               hxlabel(j) = labelitx(j,chLabels);
            end % j == 1
         elseif i == 1 && j == 2 && flgPrinting(7) == 0 % Special case row #1
            set(h,'XLim', [w_min w_max],'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'FontSize', vFontSize,'FontWeight','bold', 'FontName','Arial');
            if flgScale ==2 || flgScale == 1
               set(h,'YTick',yTick,'YTickLabel',yTickLabel, ...
                  'FontSize', vFontSize,'FontWeight','bold', ...
                  'FontName','Arial') %'Helvetica-Narrow');
            end
            if (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
               hylabel(i) = labelity(i,chLabels);
            end
         elseif j == 1 % Column #1 of lay-out
            set(h,'XLim', [w_min w_max],'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'FontSize', vFontSize,'FontWeight','bold', ...
               'FontName','Arial');
            if flgScale == 2, set(h,'YTick',yTick,'YTickLabel',yTickLabel); end
            if i == nChannels
               set(h,'FontSize', vFontSize,'FontWeight','bold', ...
                      'FontName','Arial') %'Helvetica-Narrow');
               hxlabel(j) = labelitx(j,chLabels);
            elseif flgScale == 1
             set(h,'XLim', [w_min w_max],'YLim',ylim, ...
                   'YTick',yTick,'YTickLabel',yTickLabel, ...
                   'XTick',[0 .1 .2 .3 .4 .5], ...
                   'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                   'FontSize', vFontSize,'FontWeight','bold', ...
                   'FontName','Arial');
            end
            if (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
               hylabel(i) = labelity(i,chLabels);
            end
         else  % j~=1 or j~=nChannels or not(i==1 and j==2 and ~flgPrinting(7))
            set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                  'XTick',[0 .1 .2 .3 .4 .5], ...
                  'XTickLabel',[' ';' ';' ';' ';' ';' ']);
            % Review: cases of flgScale 1 or 2 look similar.
            if flgScale == 2, set(h,'YTick',yTick,'YTickLabel',yTickLabel);end
            if flgScale == 1, set(h,'YTick',yTick,'YTickLabel',yTickLabel);end
            set(h,'FontSize', vFontSize,'FontWeight','bold', ...
                  'FontName','Arial');
         end % i == nChannels
         % Label axes; set x- y-axis limits, x- yTickLabels, FontSize

         if ~isOctave()
            h.YRuler.TickLabelGapOffset = vTickLabelGapOffset;
         end
         if nChannels > 20  % Rotate ylabels by 20 degrees if numChannels > 10.
            set(hylabel(i),'Rotation',70);
         end
         hold on

% ==============================================================================
% Plotting lower and upper bound limits of confidence interval
% ==============================================================================
         if flgPrinting(4)
            plotCI(flgPrinting,flgSignif,w,r,s,Ltmp,L2vupper,L2vlower, ...
                                          nPlotPoints,indexNotSignif, ...
                                          vPatnaikLineWidth,versionNumberMatlab)
            set(gca, ...
               'Box'         , 'on'     , ...
               'TickDir'     , 'in'     , ...
               'TickLength'  , [.02 .02] , ...
               'XMinorTick'  , 'off'      , ...
               'YMinorTick'  , 'off'      );
         end
         if flgPrinting(6) %|Coh|^2 plot in gray-lines
            plot(w,Coh2tmp,'-','LineWidth',vLineWidth,'Color',colorCoh2);
         end
         if flgPrinting(1)
            plot(h,w,Ltmp,'Color',atrib2,'LineWidth',vLineWidth);
            if c.alpha ~= 0
               plot(h,w,LTratmp,'Color',atrib1,'LineWidth',vLineWidth);
            end
         end

         if flgPrinting(2)
            atrib = 'k--'; % Patnaik significance level in black line
            plot(w,abs(getCij(Lpatnaik,r,s,nPlotPoints)),atrib, ...
               'LineWidth',vPatnaikLineWidth);
         end

% ==============================================================================
% Print GCT p-values above subplot
% ==============================================================================
%
         if r ~= s
            if (flgPrinting(1) == 1) && (flgPrinting(5) > 0)
               vaxis = double(axis);
               if isfield(c,'Tragct')
                  if Tragct(r,s)
                     if (flgPrinting(5) > 1)
                        plot(0.96*vaxis(2),0.93*(vaxis(4)-vaxis(3))+vaxis(3),'.', ...
                             'Color',vMarkerColor,'MarkerSize',vMarkerSize);
                     end
                     if isfield(c,'pvaluesgct') && (flgPrinting(5) == 1 || ...
                                                    flgPrinting(5) == 3)
                        if pvaluesgct(r,s) < 1.0e-15
                           text(vaxis(2),1.1*(vaxis(4)-vaxis(3))+vaxis(3),['<1.0e-15'], ...
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
                        text(vaxis(2),1.1*(vaxis(4)-vaxis(3))+vaxis(3), ...  % 0.6*vaxis(2)
                             sprintf('%0.3g',pvaluesgct(r,s)), ...
                             'FontSize',pValuesFontSize,'FontWeight','bold', ...
                             'Color',[0 .4 0],'HorizontalAlignment','right')
                     end
                  end
               end
            end  %strrep(num2str(pvaluesgct), '0.', '.')
         end % GCT p-values

      end % PDC2 and squared coherence
   end  % Subplot-layout row target channels loop for i = 1:nChannels
end % Subplot-layout column source channels loop for j = 1:nChannels

% ------------------------------------------------------------------------------
% Main x-axis labeling
%
supAxes = [.08 .075 .84 .84]; % [.12 .11 .84 .80];
if fs == 1,
   strXsuplabel = 'Frequency';
else
   strXsuplabel = 'Frequency (Hz)';
end

[ax1,h1] = suplabel(strXsuplabel,'x'); %,supAxes);
pos = get(ax1,'Position');

switch computer
   case 'PCWIN64'
      set(h1,'FontWeight','bold', 'FontSize',12,'FontName','Arial');
%       pos(2) = -0.03;
      pos(2) = pos(2) + 0.035; %dpos2           % 0.0545 % 0.025
   case {'GLNXA64','x86_64-pc-linux-gnu','i686-pc-linux-gnu'}
      set(h1,'FontWeight','bold', 'FontSize',12,'FontName','Arial');
%       dpos2 = 0.021*(30 - nChannels)/25; % xlabel y-shift correctionn constant
      pos(2) = pos(2) + 0.025; %dpos2           % 0.0545 % 0.025
   case 'MACI64'
      set(h1,'FontWeight','bold', 'FontSize',14,'FontName','Arial');
      %        dpos2 = 0.021*(30 - nChannels)/25; % xlabel y-shift correctionn constant
      %      pos2 = pos(2)
      pos(2) = pos(2) + 0.010;
      %       pos(2) = pos(2) + 0.025;
      %  pos(2) = pos(2) - .025; % dpos2/3;
      % save maci64_suplabelx h1 pos
   otherwise
      disp('otherwise')
      %nop
end

if ~isOctave()   
   set(ax1,'Position',pos);    % Adjust xlabel position
else
   set(ax1,'Position',pos);    % Actually suplabel does not work in Octave   
end
% ------------------------------------------------------------------------------
if flgPrinting(1) % PDC2 or DTF2

   % Main y-axis type of measure DFT, PDC, Coh or SS labeling
   %
   if strcmp(flgType,'DTF')
      gType = 'gamma';
      switch metric
         case 'euc'
            pType = 'DTF';
         case 'diag'
            pType = 'DC';
         case 'info'
            pType = '_{i}DTF';
         case 'ratio'
            pType = 'Ratio';
         otherwise
            error('Unknown metric.')
      end
   else
      gType = 'pi';
      switch metric
         case 'euc'
            pType = 'PDC';
         case 'diag'
            pType = '_{g}PDC';
         case 'info'
            pType = '_{i}PDC';
         case 'ratio'
            pType = 'iPDC/gPDC Ratio';
        otherwise
            error('Unknown metric.')
      end
   end

   switch lower(flgPrint)
      case 'screen'
         switch metric
            case 'euc'
               [ax2,h2] = suplabel(['{\mid\' gType ...
                             '_{\it{{i}{j}}}{(\lambda)\mid}^{2}}'],'y',supAxes);
            case 'diag'
               [ax2,h2] = suplabel(['{\mid{_{\it{g}}}\' gType ...
                             '_{\it{{i}{j}}}{(\lambda)\mid}^{2}}'],'y',supAxes);
            case 'info'
               [ax2,h2] = suplabel(['{\mid{_i}\' gType ...
                             '_{\it{{i}{j}}}{(\lambda)\mid}^{2}}'],'y',supAxes);
            case 'ratio'
               [ax2,h2] = suplabel(['i' flgType ' to g' flgType ' Ratio'],...
                                  'y',supAxes);
            otherwise
               error('Unknown metric.')
         end
      otherwise %Print
         [ax2,h2] = suplabel(['{|' pType '(\lambda)|}^{2}'],'y',supAxes);
         set(h2,'FontWeight','bold', 'FontSize',12);
         pos = get(ax2,'Position');
         pos(1) = pos(1) + 0.030;   % 0.0545
         set(ax2,'Position',pos);    % Adjust ylabel position
   end
elseif flgPrinting(6)
   switch lower(flgPrint)
      case 'screen'
         [ax2,h2] = suplabel(['{|\textbf{Coh}_{\it{{i}{j}}}{(\lambda)|}^{2}}'],...
                            'y',supAxes);
      otherwise %Print
         [ax2,h2] = suplabel(['{|Coh(\lambda)|}^{2}'],'y',supAxes);
         set(h2,'FontWeight','bold', 'FontSize',12);
         pos = get(ax2,'Position');
         pos(1) = pos(1) + 0.020;   % 0.0545
         set(ax2,'Position',pos);    % Adjust ylabel position
   end
elseif flgPrinting(7)
   switch lower(flgPrint)
      case 'screen'
         [ax2,h2] = suplabel(['{|SS(\lambda)|}'],'y',supAxes);
      otherwise %Print
         [ax2,h2] = suplabel(['SS(\lambda)|'],'y',supAxes);
         set(h2,'FontWeight','bold', 'FontSize',12);
         pos = get(ax2,'Position');
         pos(1) = pos(1) + 0.020;   % 0.0545
         set(ax2,'Position',pos);    % Adjust ylabel position
   end
end

%%
% Change y-LABEL font size
switch computer
   case 'PCWIN64'
      %nop
   case {'GLNXA64','x86_64-pc-linux-gnu','i686-pc-linux-gnu'}

   case 'MACI64'
      set(h2,'FontSize',16)
end


% ==============================================================================
% Fine adjustment of axis labels positions.
% ==============================================================================
%
for k = 1:nChannels
   set(hxlabel(k),'Units','normalized');

   if (flgPrinting(1) ~= 0 || flgPrinting(6) ~= 0)
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

%% labelitx
%      x-axis labeling function
%
function [hxlabel] = labelitx(j,chLabels) % Labels x-axis plottings
if isOctave()
   vFontSize = 12;
else
   switch computer
      case 'PCWIN64'
         vFontSize = 10;
      case {'GLNXA64','x86_64-pc-linux-gnu','i686-pc-linux-gnu'}
         vFontSize = 10;
      case 'MACI64'
         vFontSize = 12;
   end
end

if isempty(chLabels)
   hxlabel = xlabel(['j=' int2str(j)]);
   set(hxlabel,'FontSize',vFontSize, 'FontWeight','bold', ...
               'FontName','Arial') % 'FontName','Helvetica-Narrow'
else
   hxlabel = xlabel([chLabels{j}]);
   set(hxlabel,'FontSize',vFontSize,'FontWeight','bold', ...
               'FontName','Arial');
end

pos = get(hxlabel,'Position');
pos(2) = pos(2) + 0.085;   % 0.0545
set(hxlabel,'Position',pos);

end

%% latelity 
% y-axis labeling function
%
function [hylabel] = labelity(i,chLabels) % Labels y-axis plottings
if isOctave()
   vFontSize = 12;
else
   switch computer
      case 'PCWIN64'
         vFontSize = 10;
      case {'GLNXA64','x86_64-pc-linux-gnu','i686-pc-linux-gnu'}
         vFontSize = 10;
      case 'MACI64'
         vFontSize = 12;
   end
end
if isempty(chLabels)
   hylabel = ylabel(['i=' int2str(i)],...
      'Rotation',90);
   set(hylabel,'FontSize',vFontSize,'FontWeight','bold', ...
               'FontName','Arial') %      'FontName','Helvetica-Narrow')  % 'Times'
else
   hylabel = ylabel([chLabels{i}]);
   set(hylabel,'FontSize',vFontSize,'FontWeight','bold');
end

pos = get(hylabel,'Position');
if isOctave()
   pos(1) = pos(1) - 0.0545;   % 0.0545
else
   pos(1) = pos(1) + 0.0845;   % 0.0545
end
set(hylabel,'Position',pos)
end

%% plotci
%      Plot upper and lower confidence interval limits
%
function plotCI(flgPrinting,flgSignif,w,r,s,Ltmp,L2vupper,L2vlower, ...
               nPlotPoints,indexNotSignif,vPatnaikLineWidth,versionNumberMatlab)

if flgPrinting(4) >= 2 % Plot error bar
   if flgPrinting(4) == 3
      kstep = 2; % Errorbar spacing parameter
   elseif flgPrinting(4) == 4
      kstep = round(length(w)/32); % Limiting number of error bar on x-axis
      %                     % to 32.
   elseif flgPrinting(4) == 5
      kstep = round(length(w)/16); % Limiting number of error bar on x-axis
      %                     % to 16.
   end

   L2vuppertmp = abs(getCij(L2vupper,r,s,nPlotPoints));
   L2vuppertmp(indexNotSignif) = NaN;
   L2vlowertmp = abs(getCij(L2vlower,r,s,nPlotPoints));
   L2vlowertmp(indexNotSignif) = NaN;

   Ltmp2 = Ltmp; %= = abs(getCij(L,r,s,nPlotPoints));
   Ltmp2(indexNotSignif) = NaN;
   varPDC2upper = L2vuppertmp - Ltmp2;

   % Sparsed error bar plotting for confidence interval
   if flgPrinting(4) == 2
      web = 1:length(w);
   else
      web = kstep:kstep:length(w); % for sparse plotting interval.
   end
   b = zeros(length(web),2,1);
   b(:,1,1) = varPDC2upper(web);
   b(:,2,1) = varPDC2upper(web);

   Nweb = length(web);

   if flgPrinting(4) == 2  % Confidence interval with shadeplotting
      web = 1:length(w);
      indexTemp = ones(1,Nweb);
      %web(webb) index with NaN values.
      isnanIndex = find(isnan(varPDC2upper(web)));
      indexTemp(isnanIndex) = 0;
      webb = web(indexTemp > 0);      % webb contains not NaN web indices.
      webbProbe = [-1 webb Nweb+2]; % Temporary variable for probing not NaN elements.
      webbDiff = diff(webbProbe);
      indexNotNaN = find(webbDiff > 1);
      %webbDiffIndex = webb(webbDiff(1:end-1) > 1);
      nDiff2 = sum(webbDiff > 1);

      if nDiff2 ~= 0
         for kDiff = 1:nDiff2-1
            indexPlotVar0   = webb(indexNotNaN(kDiff));
            indexPlotVarEnd = webb(indexNotNaN(kDiff+1)-1);

            webbIndexTemp = web(indexPlotVar0:indexPlotVarEnd);
            wtmp = w(webbIndexTemp);     wtmp = [wtmp(1) wtmp];
            ytmp = Ltmp(webbIndexTemp);  ytmp = [ytmp(1); ytmp];
            btmp = b(webbIndexTemp,:,:); btmp = [[0 0]; btmp];
            hE = boundedline(wtmp,ytmp,btmp,'-k','alpha'); % Transparent patch
         end
      else
         wtmp = w(webb);     wtmp = [wtmp(1) wtmp];
         ytmp = Ltmp(webb);  ytmp = [ytmp(1); ytmp];
         btmp = b(webb,:,:); btmp = [[0 0]; btmp];
         hE = boundedline(wtmp,ytmp,btmp(:,1),'-k','alpha'); % Transparent patch
      end

   else % Confidence interval with error bars
      if isOctave()
         hE = errorbar(w(web),Ltmp(web),varPDC2upper(web));
      else
         if versionNumberMatlab > 9.0 % > 2016a: 'CapSize' appeared on 2016b.
            hE = errorbar(w(web),Ltmp(web),varPDC2upper(web), ...
                                                 'LineWidth',0.75,'CapSize',4);
         else % for older v. Matlab and Octave
            hE = errorbar(w(web),Ltmp(web),varPDC2upper(web), ...
                                                             'LineWidth',0.75);
         end
      end
   end

   if exist('hE','var')
      set(hE,'LineWidth',1,'color',[0 0 0]); %[0.4 0.4 0.4]);
   end

else  % Plot confidence interval with gray dashed-line
   if flgPrinting(4) && flgSignif
      atrib = 'k--';

      % Lower-bound confidence interval limits
      tmp = getCij(L2vlower,r,s,nPlotPoints);
      % lower-bound can be negative.
      tmp(indexNotSignif) = NaN;
      plot(w,tmp, atrib,'LineWidth',vPatnaikLineWidth,'Color',[.4 .4 .4]);

      % Upper-bound confidence interval limits
      tmp = abs(getCij(L2vupper,r,s,nPlotPoints));
      tmp(indexNotSignif) = NaN;
      plot(w,tmp, atrib,'LineWidth',vPatnaikLineWidth,'Color',[.4 .4 .4]);
   end
end
end  % Plotting lower and upper bound of confidence interval
