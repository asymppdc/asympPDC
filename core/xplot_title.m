%% XPLOT_TITLE
%        Put a title above PDC/DTF matrix layout subplots with measure and
%        alpha-value information.
%
%% Syntax
%   XPLOT_TITLE(c)
%   XPLOT_TITLE(alpha,metric,measure)
%
%% Input Arguments
%   alpha:      significance level
%
%   measure:    should be either 'pdc'| 'dtf'
%
%   metric:     'euc'  -- Euclidean   ==> |PDC|^2  or |DTF|^2
%               'diag' -- diagonal    ==> |gPDC|^2 or |DC|^2
%               'info' -- information ==> |iPDC|^2 or |iDTF|^2
%
%% Output 
%   Print a title above all subplots of DTF/PDC.
%
% See also XPLOT, PVALUES_XPLOT
%          xplot | <xplot.html> |

% (C) Koichi Sameshima & Luiz A. Baccala, 2021. See file license.txt in
% installation directory for licensing terms.
%

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

if isOctave()
   set(h,'FontSize',14)   
else
   set(h,'FontSize',10)
end