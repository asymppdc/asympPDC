function [Tr_igct, pValue_igct] = igct_alg(u,A,pf,igct_signif,flgVerbose)
%IGCT_ALG  Perform instantaneous Granger causality test (IGCT) from time series,
%          AR model coefficients and covariance matrix.
%
% Syntax:
%        [Tr_igct, pValue_igct] = IGCT_ALG(u,A,pf,igct_signif,flgVerbose)
%
% Input arguments:
%        u           - multivariate time series data
%        A           - Autoregressive model coefficients matrix
%        pf          - Covariance matrix
%        igct_signif - Significance level for null hypothesis for IGCT
%        flgVerbose  - If 1: verbose for printing the testing results on screen
%                        0: perform tests silently
%
% Output arguments:
%        Tr_igct     - IGCT based connectivity matrix
%        pValue_igct - IGCT p-values matrix
%
% Reference:
% [1] Lutkepohl, H (2005). New Introduction to Multiple Time Series Analysis. 
%                         Springer-Verlag. 
%
% See also: GCT_ALG 

% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.


   if nargin < 4
      error('IGCT_ALG.M requires at least 4 input arguments.');
   elseif nargin == 4
      flgVerbose = 0;
   end

   % Checking input parameters
   [nChannels,~,IP] = size(A);
   nSegLength = length(u);
   Z = zmatrm(u,IP); 
   gamma = Z*Z';

   idx = eye(nChannels)==1;

   % Instantaneous Granger causality test routine
   %
   [Tr_igct,Va_igct,v_igct,th_igct,pValue_igct] = granmaty(pf,nSegLength, ...
                                                                  (1-igct_signif));

   %  Main diagonal elements are filled with NaN.
   Tr_igct(idx) = NaN;
   pValue_igct(idx) = NaN;

   if flgVerbose
      disp(repmat('=',1,100))
      disp('                  INSTANTANEOUS GRANGER CAUSALITY TEST')
      fprintf(repmat('-',1,100))
      fprintf('\n')
      disp('Instantaneous connectivity matrix:')
      disp(Tr_igct)
      fprintf('\n')
      fprintf('Instantaneous Granger Causality test p-values:')
      fprintf('\n')
      disp(pValue_igct)
      nPairsIGC = (sum(sum(Tr_igct==1)))/2;

      if nPairsIGC == 0
         fprintf('\n>>>> Instantaneous Granger Causality NOT detected.\n')
      elseif nPairsIGC == 1
         fprintf('\n>>>> There is a pair of channels with significant Instantaneous ')
         fprintf('\n     Granger Causality.\n')
      else
         fprintf(['\n>>>> There are ' int2str(nPairsIGC) ' pairs of channels with'])
         fprintf('\n      significant Instantaneous Granger Causality.\n')
      end
      fprintf('\n')
      disp(repmat('=',1,100))
   end
   end

% Subfunctions:

   %==========================================================================
   function [Tr,Va,v,th,pValue]=granmaty(pf,N,significance)
   % Test Granger causality structure
   %
   %[Tr,Va,v,th,pValue]=granmaty(SU,N,significance);
   % Program to test granger causality structure
   %
   % input: N (number of points)
   %        pf - covariance of modelling errors
   %        significance - test significance level
   %
   % output: Tr -test result matrix (i,j) entry=1 j->i causality cannot
   %             be rejected
   %         Va - test value matrix
   %         v  - degrees of freedom
   %         th - threshold value
   %         pValue - test p-value
   %
   % % 01/30/1998 - L.A.B.
   % % 27/10/2009 - Stein - Mudou para v=1.
   %
   % disp('Instantaneous Granger causality test: ');
   % significance
   [n, m]=size(pf);
   Va=zeros(n,m);
   Tr=zeros(n,m);
   CO=zeros(n,m);
   pValue=zeros(n,m);
   for i=1:n
      for j=1:n
         if i>j
            CO(i,j)=1;
            [Tr(i,j),Va(i,j),v,th,pValue(i,j)]=instata(CO,pf,N,significance);
            Tr(j,i)=Tr(i,j);
            Va(j,i)=Va(i,j);
            CO(i,j)=0;
            pValue(j,i)=pValue(i,j);
         end
      end
   end
   end

   %==========================================================================
   function [y,value,v,th,pValue]=instata(CO,pf,N,significance)
   % Test for instataneous causality
   % input: CO - matrix describing the structure for testing - 1 position to test.
   %        pf - residual covariance
   %        N - number of poinst
   %
   % output: y - test result - 0 instantaneous causality rejected - 1 not rejected
   %         value - test value
   %         v - degrees of freedom # constraints.
   %         th -threschold
   si=vech(pf);
   CO=tril(CO);
   [m,n]=size(CO);
   lb=length(si);
   Ct=vech(CO);
   Ct1=zeros(size(Ct'));
   Ctf=[ ];
   l=sum(Ct');
   for i=1:length(Ct)
      if Ct(i) == 1
         Ct1(i) = 1;
         Ctf = [Ctf; Ct1];
         Ct1 = zeros(size(Ct'));
      end
   end
   C=Ctf;
   ln = length(pf);
   D = pinv(dmatrix(ln));
   value = N*(C*si)'*inv(2*C*D*kron(pf,pf)*D'*C')*C*si;
   v = 1; %2; Chi-square distribution degree of freedom. C.S. Changed to v = 1.
   th = chi2inv(significance,v);
   y = value >= th;
   pValue = 1-chi2cdf(value,v); % p-value of instantaneous Granger causality test
   end


   %==========================================================================
   %
   %  01/30/1998 - L.A.B.
   %
   function D=dmatrix(m)
   D=zeros(m*m,m*(m+1)/2);
   u=[ ];
   v=[];
   for j=1:m
      for i=1:m
         u=[u ;[i j]];
         if j<=i
            v=[v ;[i j]];
         end
      end
   end
   w=fliplr(v);
   for i=1:m*m
      for j=1:m*(m+1)/2
         if sum(u(i,:)==v(j,:))==2
            D(i,j)=1;
         end
      end
      for j=1:m*(m+1)/2
         if sum(u(i,:)==w(j,:))==2
            D(i,j)=1;
         end
      end
   end
   end


   %==========================================================================
   % VECH or VEC is matrix column stacking operator function
   %
   %function y=vech(Y);
   %
   % input:  Y - matrix
   % output: y - Stacked column vector
   %
   % % 01/30/1998 - L.A.B.

   function y=vech(Y)
   y = [ ];
   [m, n] = size(Y);
   for i = 1:m
      y = [y ;Y(i:n,i)];
   end
   end


   %==========================================================================
   % Computation of Z - data structure (no estimation of the mean)
   %
   % function Z=zmatr(Y,p);
   %
   % input:  Y - data in row vectors
   %         p - model covariance order
   %
   % output: Z
   %
   % 01/30/1998 - L.A.B.
   %
   function Z=zmatrm(Y,p)
   [K, T] = size(Y);
   y1 = [zeros(K*p,1);reshape(flipud(Y),K*T,1)];
   Z =  zeros(K*p,T);
   for i = 0:T-1
      Z(:,i+1) = flipud(y1(1+K*i:K*i+K*p));
   end
   %Z=[ones(1,T);Z];
   end
