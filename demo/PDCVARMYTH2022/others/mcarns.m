%% MCARNS 
%       Nuttall-Strand algorithm for autoregressive model estimation.
%
%% Syntax:
%       [pf,A,pb,B,ef,eb] = MCARNS(u,IP)
%
%% Input Arguments
%     IP     - Order of autoregressive model (integer)
%     u      - Complex matrix with NUMCHS channels of sample data
%
%%  Output Arguments:
%     pf     - Covariance matrix of NUMCHS x NUMCHS of linear forward
%              prediction error
%     A      - Complex array of forward linear prediction matrix
%              coefficients
%     pb     - Complex backward linear prediction error covariance array
%     B      - Complex array of backward linear prediction matrix
%              coefficients
%     ef     - Forward residuals
%     eb     - Backward residuals
%
%% Description:
%   This function calulate the coeficients of multi-channel autoregressive 
%   matrix using Nuttall-Strand algorithm (a generalization of single channel 
%   harmonic method).
%
%% Notes
%   This MATLAB implementation of MCARNS is a translation of FORTRAN code from 
%   Appendix 15.B page 424 in [1]. (KS, 1998)
%
%% References:
%   [1] Marple Jr SL (1987). Digital Spectral Analysis with Application.
%       Prentice-Hall, Englewood-Cliffs. 
% 
%  Equation numbering refers to [1].
%

function [pf,A,pb,B,ef,eb,ISTAT]=mcarns(u,IP)

if nargin ~= 2, error('MCARNS requires 2 input parameters.'); end

[lx,cx] = size(u);

if lx > cx, error('Input matrix is probably transposed.'), end
NUMCHS = lx;          % Number of channels.
MAXORDER = 200;       % Maximum order of AR model allowed for calculation.
N = max(size(u));     % N - Number of samples per channel.

%   ** Initialization **
if (IP > MAXORDER)
   error('IP > 200.');
end

ef = u;        % Eq. (15.91)
eb = u;        % Eq. (15.91)

pf = u*u';     % Eq. (15.90)
pb = pf;       % Eq. (15.90)

M = 0;
%   ** Main Loop **
while 1
   %  Update estimated covariance errors               Eq. (15.89)
   pfhat = ef(:,M+2:N)*ef(:,M+2:N)';
   pbhat = eb(:,M+1:N-1)*eb(:,M+1:N-1)';
   pfbhat = ef(:,M+2:N)*eb(:,M+1:N-1)';

   M = M+1;

   %  Calculate estimated partial correlation matrix - Eq. (15.98)
   %             (Nuttall-Strand algorithm only)
   RHO = lyap(pfhat*inv(pf),inv(pb)*pbhat,-2*pfbhat);

   %  Update forward and backward reflection coeficients
   %  Eqs. (15.73),(15.74),(15.78) (algoritjm  by Nuttall-Strand)
   AM = -RHO*inv(pb);
   BM = -RHO'*inv(pf);

   A(:,:,M) = AM;
   B(:,:,M) = BM;
   %
   %  Update forward and backward covariance error  - Eqs. (15.75),(15.76)
   pf = pf-AM*BM*pf;
   pb = pb-BM*AM*pb;

   %  Update forward and backward predictor coeficients - Eqs.(15.84),(15.85)
   if M ~= 1
      for K = 1:M-1
         temp1 = A(:,:,K);
         A(:,:,K) = A(:,:,K)+AM*B(:,:,M-K);
         B(:,:,M-K) = B(:,:,M-K)+BM*temp1;
      end
   end
   %  Update residuals
   Tef = ef;
   ef(1:NUMCHS,N:-1:M+1) = ef(:,N:-1:M+1)+AM*eb(:,N-1:-1:M);
   eb(1:NUMCHS,N:-1:M+1) = eb(:,N-1:-1:M)+BM*Tef(:,N:-1:M+1);

   %  Verify if model order is adequate
   if M == IP, A = -A; B = -B; break, end
end
