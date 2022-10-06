%% MVARRESIDUE
%        Portmanteau residues test for whiteness
%
%% Syntax:
%        [Pass,Portmanteau,st,ths] = MVARRESIDUE(ef,ns,p,aValue,h,flgVerbose)
%
%% Input Arguments:
%        ef       - matrix of column vectors - residues
%        ns       - number of points estimated
%        p        - model order
%        aValue   - confidence interval
%        h        - lag  - maximum lag
%        flgVerbose - If 1: verbose for printing the testing results
%                        0: perform tests silently
%
%% Output Arguments: 
%         Pass    - number of computed correlation coefficients
%         Portmanteau -result of Portmanteau test: 0 reject; 1 not rejected
%                                                           (white hypothesis)
%         st      - Portmanteau statistic
%         ths     - threshold value
%
%% Reference:
% [1] Lutkepohl, H (2005). New Introduction to Multiple Time Series Analysis.
%                         Springer-Verlag. Sec. 4.4.3, p. 169--171
%
%% See also: mvar
%
%         LAB 11/Apr/2000
%         KS  28/Apr/2007

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.

%%

function [Pass,Portmanteau,st,ths]=mvarresidue(ef,ns,p,aValue,h,flgVerbose)

if ~exist('h','var'), h=10; end
ef=ef';
[t,Portmanteau,st,ths,compr]=crosstest(ef,h,ns,aValue,p);
Pass=t/compr;

if flgVerbose
   fprintf(['\n' repmat('=',1,100) '\n'])
   disp('                  MVAR RESIDURES TEST FOR WHITENESS')
   disp(repmat('-',1,100))
end

if flgVerbose
   if Portmanteau
      disp(['Good MAR model fitting! Residues white noise hypothesis ' ...
         'NOT rejected.'])
   else
      disp('(**) Poor MAR model fitting:')
      disp('                   Residues white noise hypothesis was rejected.')
   end
   disp(['Pass = ', sprintf('%g',Pass)])
   disp(['  st = ', sprintf('%g',st)])
end
end


%% crosstest
%      Cross-correlation computation between data signal column vectors
%
%% Syntax
%          [t,Portmanteau,s,ths,compr,X] = crosstest(u,h,ns,th,p);
%
%% Input arguments
%        u = matrix of column vectors - residues
%        h - lag - maximu lag
%        ns - number of points estimated
%        th - confidence
%        p - model order
%
%% Output arguments
%         t - correlation test - number of times the 2/sqrt(ns) threshold is exceeded
%         Portmanteau - test result - 0 - REJECT , 1 not rejected (white hypothesis)
%         s - portamanteau statistic
%         ths - threshold value
%         compr - # of computed correlation coefficients
%         X - cross correlation vector columns: order 11 12 ... 21 22 .. etc
%                   (i-1)*nseries+j
%         reference:Luktpohl(1993) Chap. 4
%
%         LAB 11/4/2000
%         Stein 27/10/2009 - Changes (h -> h+1 and ns-i -> ns-i+1)
%

function [t,Portmanteau,s,ths,compr,X] = crosstest(u,h,ns,th,p)
if nargin==4
   p=0;
end
[~,n]=size(u);
C=zeros(n,n,h);
X=xcorr(u,h,'coeff');
for i=h+1:2*h+1
   C(:,:,i-h)=reshape(X(i,:),n,n);
end
%Simple test
t=sum(sum(sum(abs(C(:,:,2:h+1))>2/sqrt(ns))));
% Portmanteau
s=0;
SO=pinv(C(:,:,1));
for i=2:(h+1) %Stein - changed to h+1
   s=s+((ns-i+1)^(-1))*trace(C(:,:,i)'*SO*C(:,:,i)*SO); %Stein - changed to ns-i+1
end
s=ns^2*s;
ths=chi2inv(th,(h-p)*n^2);
Portmanteau=s<ths;
compr=n*n*h;
end