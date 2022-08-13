%% GCT_ALG
%        Perform Granger causality test (GCT) from time series, AR model
%        coefficients and covariance matrix.
%
%% Syntax:
%        [Tr_gct, pValue_gct] = GCT_ALG(u,A,pf,gct_signif,igct_signif,
%                                                                    flgVerbose)
%
%% Input arguments:
%        u           - data
%        A           - Autoregressive model coefficient matrix
%        pf          - Covariance matrix
%        gct_signif  - Significance level for null hypothesis for GCT
%        igct_signif - Significance level for null hypothesis for IGCT
%        flgVerbose  - If 1: verbose for printing the testing results
%                         0: perform tests silently
%
%% Output arguments:
%        Tr_gct      - GCT based connectivity matrix
%        pValue_gct  - GCT p-values matrix
%
%% Reference:
% [1] Lutkepohl, H (2005). New Introduction to Multiple Time Series Analysis. 
%                         Springer-Verlag. 
%
%% See also: IGCT_ALG

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.


function [Tr_gct, pValue_gct] = gct_alg(u,A,pf,gct_signif,flgVerbose)
%%
% Checking input parameters

if nargin < 4
   error('GCT_ALG.M requires at least 3 input parameters.');
elseif nargin == 4
   flgVerbose = 0;
end


%%
%

[nChannels,~,IP] = size(A);
nSegLength = length(u);
Z = zmatrm(u,IP); 
gamma = Z*Z';

idx = eye(nChannels)==1;

%%
% Granger causality test routine
b = reshape(A,nChannels*nChannels*IP,1);
[Tr_gct,Va_gct,v_gct,th_gct,pValue_gct] = granmatx(b,gamma,pf,(1-gct_signif));

%%
% Main diagonal elements are filled with NaN.
Tr_gct(idx) = NaN;
pValue_gct(idx) = NaN;

if flgVerbose,
   format compact
   %   disp('Granger causality test:')
   fprintf('\n')
   disp(repmat('=',1,100))
   disp('                         GRANGER CAUSALITY TEST')
   fprintf(repmat('-',1,100))
   fprintf('\n')
   disp('Connectivity matrix:')
   disp(Tr_gct)
   fprintf('\n')
   disp('Granger causality test p-values:')
   disp(pValue_gct)
   fprintf('\n')
   %    disp('Causality value:')
   %    Va_gct
end
end

function [Tr,Va,v,th,pValue]=granmatx(b,G,SU,significance)

% function [Tr,Va,v,th,pValue]=granmatx(b,G,SU);
% Routine to test Granger causality structure
%
% input: b - reshaped A matrix (A,1,Numch times p)
%        G - data covariance of vector times N (number of points)
%        SU - covariance of modelling errors
%        significance - statistical significance level
%
% output: Tr - test result matrix (i,j) entry=1 j->i causality cannot
%              be rejected
%         Va - test value matrix
%         v  - degrees of freedom
%         th - threshold value
%         pValue - p-value
% % 01/30/1998 - L.A.B.
%

[n,m] = size(SU);
Va = zeros(n,m);
Tr = zeros(n,m);
CO = zeros(n,m);
pValue = zeros(n,m);

for i=1:n
   for j=1:n
      if i~=j
         CO(i,j)=1;
         [Tr(i,j),Va(i,j),v,th,pValue(i,j)]=grangt(CO,b,G,SU,significance);
         CO(i,j)=0;
      end
   end
end
end


% =========================================================================

function [y,value,v,th,pValue]=grangt(CO,b,G,SU,significance)
%%% Function GRANGT for Granger Causality Test
%function [y,value,v,th,pValue]=grangt(CO,b,G,SU,significance);
%
% Causality test
%
% input: CO - matrix describing the structure for testing - 1 position to test.
%        b - parameter vector
%        G - Gamma*T - data covariance matriz times T record length
%        SU - residual covariance
%        significance - statistical significance level
%
% output: y - test result - 0 granger causality rejected - 1 not rejected
%         value - test value
%         v - degrees o freedom # oconstraints.
%         th -threschold
%
% Test for GrangerCausality
%
% % 01/30/1998 - L.A.B.
%

[n,m] = size(CO);
lb = length(b);
p0 = lb/(m*n);
Ct = reshape(CO,1,m*n);
Ct1 = [];
for i = 1:p0
   Ct1 = [Ct1 Ct];
end
Ct = Ct1;
l = sum(Ct);
K = zeros(l,m*n*p0);
for i = 1:l
   [p ,q] = max(Ct);
   K(i,q) = 1;
   Ct(q) = 0;
end
C = K;

value = (C*b)'*inv(C*kron(inv(G),SU)*C')*C*b;
v = l;
th = chi2inv(significance,v);

y = value >= th;
pValue = 1 - chi2cdf(value,v);
end


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

[n,m]=size(pf);
Va=zeros(n,m);
Tr=zeros(n,m);
CO=zeros(n,m);
pValue=zeros(n,m);
for i=1:n
   for j=1:n
      if i>j
         CO(i,j)=1;
         [Tr(i,j),Va(i,j),v,th,pValue(i,j)] = instata(CO,pf,N,significance);
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

si = vech(pf);
CO = tril(CO);
[m,n] = size(CO);
lb = length(si);
Ct = vech(CO);
Ct1 = zeros(size(Ct'));
Ctf = [ ];
l = sum(Ct');
for i = 1:length(Ct)
   if Ct(i) == 1
      Ct1(i) = 1;
      Ctf = [Ctf; Ct1];
      Ct1 = zeros(size(Ct'));
   end
end
C = Ctf;
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
function D = dmatrix(m)
D = zeros(m*m,m*(m+1)/2);
u = [ ];
v = [];
for j = 1:m
   for i = 1:m
      u = [u ;[i j]];
      if j <= i
         v = [v; [i j]];
      end
   end
end
w = fliplr(v);
for i = 1:m*m
   for j = 1:m*(m+1)/2
      if sum(u(i,:) == v(j,:))==2
         D(i,j) = 1;
      end
   end
   for j = 1:m*(m+1)/2
      if sum(u(i,:) == w(j,:))==2
         D(i,j) = 1;
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

function y = vech(Y)
y = [ ];
[m,n] = size(Y);
for i = 1:m
   y = [y; Y(i:n,i)];
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
[K,T] = size(Y);
y1 = [zeros(K*p,1);reshape(flipud(Y),K*T,1)];
Z =  zeros(K*p,T);
for i = 0:T-1
   Z(:,i+1) = flipud(y1(1+K*i:K*i+K*p));
end
end
