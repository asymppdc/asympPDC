%% STANDPLOTX2
%      Standard plot of matricial layout of connectivity measures
%
%% Syntax
%         STANDPLOTX2(L,w,limits,holdFlagToggle,atrib,yAxisType,aLineWidth,pColor4,aCoord)
%
%% Input arguments
%       L             - Connectivity measures to plot
%       w             - Frequency range
%       limits        - x and y-axis limits
%       flgHoldToggle - hold on/off flag choice
%       atrib         - Line color
%       yAxisType     -
%       aLineWidth    - Specify line width 
%       pColor4       - NOT USED
%       aCoord        - Choose subplot coordinate to be plotted
%
%% Example
%   standplotx2(pdc,[],[0 .5 -0.75 1.25],flghold, C(4,:),lWidth(k))
%               L    w    axis limits   flgHoldToggle atrib
%

function standplotx2(L,w,limits,flgHoldToggle,atrib,yAxisType,aLineWidth)

[nChannels,nChannels,nFreqs] = size(L);
k = 0;
wx = (0:nFreqs-1)*.5/nFreqs;

atrib2 = [0.41,0.41,0.41];
chLabels = [];

if nargin < 6
   yAxisType = 0;
end
if nargin < 7
   aLineWidth = 3.0;
end
if nargin < 8
   pColor4 = 1;
end
if nargin < 9
   aCoord = ones(nChannels,nChannels);
end
if nargin<2
   w = [];
   flgHoldToggle = 0;
end
if isempty(w)
   w = wx;
end
if nargin < 5
   atrib = atrib2;
end
if nargin < 3
   flgHoldToggle = 0;
end

for i = 1:nChannels
   for j = 1:nChannels
      k = k+1;
      if aCoord(i,j)
         h = subplot2(nChannels,nChannels,k);
         y = reshape(L(i,j,1:nFreqs),nFreqs,1);
         p = plot(w, y, 'Color', atrib, 'LineWidth', aLineWidth);
         
         if yAxisType ~= 0
            axis(limits)
         end
         
         if i == nChannels
            set(h,'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' 0';'  ';'  ';'  ';'  '; '.5'], ...
               'FontSize',10)
            if j == 1
               hylabel(i) = labelity(i,chLabels);
               hxlabel(j) = labelitx(j,chLabels);
            else
               hxlabel(j) = labelitx(j,chLabels);
            end
         elseif j == 1
            hylabel(i) = labelity(i,chLabels);
            set(h, 'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'FontSize',10)
         else
            set(h, 'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'FontSize',10)
         end
         
         grid on
         if flgHoldToggle
            hold on,
         end
      end
   end
end

%==========================================================================
function [hxlabel] = labelitx(j,chLabels) % Labels x-axis plottings
if isempty(chLabels)
   hxlabel = xlabel(['j = ' int2str(j)]);
   set(hxlabel,'FontSize',12);
else
   hxlabel = xlabel([chLabels{j}]);
   set(hxlabel,'FontSize',12)
end

% ========================================================================

function [hylabel] = labelity(i,chLabels) % Labels y-axis plottings
if isempty(chLabels)
   hylabel = ylabel(['i = ' int2str(i)],...
      'Rotation',90);
   set(hylabel,'FontSize',12);
else
   hylabel = ylabel([chLabels{i}]);
   set(hylabel,'FontSize',12);
end

%==========================================================================

