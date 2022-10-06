%% XPLOT_TITLE
%        Put a title above PDC/DTF matrix layout subplots with measure and
%        alpha-value information.
%
%% Syntax:
%   XPLOT_TITLE(c)
%   XPLOT_TITLE(alpha,metric,measure)
%
%% Input Arguments:
%   alpha:      significance level
%   measure:    should be either 'pdc'| 'dtf'
%   metric:     'euc'  -- Euclidean   ==> |PDC|^2  or |DTF|^2
%               'diag' -- diagonal    ==> |gPDC|^2 or |DC|^2
%               'info' -- information ==> |iPDC|^2 or |iDTF|^2
%
%% Output: 
%   Print main title above all subplots set after xplot/xplot_pvalues execution.
%
%% See also: XPLOT, XPLOT_PVALUES

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.

function [] = xplot_title(alpha,metric,measure,obs)

if nargin == 1
   c=alpha;
   if isfield(c,'dtf2')
      measure = 'DTF';
      alpha = c.alpha;
      metric =  c.metric;
   elseif isfield(c,'pdc2')
      measure = 'PDC';
      alpha = c.alpha;
      metric =  c.metric;
   end
elseif nargin < 3
   measure = 'PDC';
   obs = '';
elseif nargin < 4
   if ~(strcmp(upper(measure),'PDC') || strcmp(upper(measure),'DTF'))
      obs = measure;
      measure = 'pdc';
   else
      obs = '';
   end
end

if isempty(measure)
   measure = 'PDC';
end

if exist('c') ~= 0
   alphastr = sprintf('%0.3g',100*(c.alpha));
else
   alphastr = sprintf('%0.3g',100*alpha);
end

if ~(strcmp(upper(measure),'PDC') || strcmp(upper(measure),'DTF'))
   strMeasure = measure;
end

switch lower(metric)
   case 'euc'
      switch upper(measure)
         case {'DTF','PDC'}
            strMeasure = ['{|' upper(measure) '|}^2'];
         otherwise
            %nop
      end
   case 'diag'
      switch upper(measure)
         case 'DTF'
            strMeasure = '{|DC|}^2';    % generalized DTF or Directed Coherence (DC)
         case 'PDC'
            strMeasure = '{|_gPDC|}^2'; % generalized PDC
         otherwise
            %nop
      end
   case 'info'
      switch upper(measure)
         case {'DTF','PDC'}
            strMeasure = ['{|_i' upper(measure) '|}^2'];
         otherwise
            %nop
      end      
   otherwise
      error('Unknown metric.')
end

if ~isempty(obs)
   obs = [' ' obs];
end

if alpha <= 0
   [ax,h] = suplabel([strMeasure obs],'t');
else
   [ax,h] = suplabel([strMeasure ' (' '{\alpha = ' alphastr '%}' ')' ...
                      obs], 't');
end

set(ax,'Units','normalized')
pos = get(ax,'Position');

if isOctave()
   set(h,'FontSize',10)
else
   
   switch computer
      case 'PCWIN64'
         set(h,'FontSize',8)
         pos(4) = pos(4) + 0.0125;
      case {'GLNXA64','x86_64-pc-linux-gnu'}
         set(h,'FontSize',8)
         pos(4) = pos(4) + 0.0100;
      case 'MACI64'
         set(h,'FontSize',10)
         pos(4) = pos(4) + 0.0080;
      otherwise
         set(h,'FontSize',8)
   end
   set(ax,'Position',pos);
end

drawnow

end
% ax : Axes (suplabel) with properties:
%              XLim: [0 1]
%              YLim: [0 1]
%            XScale: 'linear'
%            YScale: 'linear'
%     GridLineStyle: '-'
%          Position: [0.0900 0.0700 0.8550 0.8950]
%             Units: 'normalized'