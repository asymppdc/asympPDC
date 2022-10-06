%% MCARNS 
%    Nuttall-Strand algorithm for autoregressive model estimation 
%
%% Syntax:
%   [pf,A,pb,B,ef,eb] = MCARNS(u,IP)
%
%% Input Arguments:
%     u      - Complex matrix with NUMCHS channels of sample data
%     IP     - Order of autoregressive model (integer)
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
%   This function calculate the coefficients of multi-channel autoregressive
%   matrix using Nuttall-Strand algorithm (a generalization of single channel
%   harmonic method).
%
%% Notes:
%      * This implementation is a FORTRAN code translation from
%        Appendix 15.B page 424 of Marple Jr.(1987). (KS 1998)
%      * Equation numbers are identifical in references [1] and [2].
%
%% References:
%   [1] Marple Jr., SL. Digital Spectral Analysis with Application.
%       Prentice-Hall, Englewood-Cliffs, 1987. 
%   [2] Marple Jr., SL. Digital Spectral Analysis with Application. 2nd Ed.
%       Dover Publications, Inc., Mineola New York, 2019. 
%
%% See also: MVAR, MCARvm, ARTFITCAPS, CMLSM 

% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.

%%

function [pf,A,pb,B,ef,eb,ISTAT] = mcarns(u,IP)

if nargin ~= 2, error('MCARNS requires 2 input arguments.'); end

[lx,cx] = size(u);
if lx > cx, error('Input matrix is possibly transposed.'), end

NUMCHS   = lx;        % Number of channels.
MAXORDER = 200;       % Maximum order of AR model allowed for calculation.
N = max(size(u));     % N - Number of samples per channel.

%   ** Initialization **
if (IP > MAXORDER)
   error('IP > 200.');
end

ef = u;        % Eq. (15.91)
eb = u;        % Eq. (15.91)

pf = u * u';   % Eq. (15.90)
pb = pf;       % Eq. (15.90)

m = 0;

%   ** Main Loop **
while 1
   %  Update estimated covariance errors               Eq. (15.89)
   pfhat  = ef(:,m+2:N)   * ef(:,m+2:N)';
   pbhat  = eb(:,m+1:N-1) * eb(:,m+1:N-1)';
   pfbhat = ef(:,m+2:N)   * eb(:,m+1:N-1)';

   m = m + 1;

   %  Calculate estimated partial correlation matrix - Eq. (15.98)
   %             (Nuttall-Strand algorithm only)
   RHO = lyap(pfhat*inv(pf),inv(pb)*pbhat,-2*pfbhat);

   %  Update forward and backward reflection coefficients
   %  Eqs. (15.73),(15.74),(15.78) (algorithm  by Nuttall-Strand)
   AM = -RHO  * inv(pb);
   BM = -RHO' * inv(pf);

   A(:,:,m) = AM;
   B(:,:,m) = BM;
   
   %  Update forward and backward covariance error  - Eqs. (15.75),(15.76)
   pf = pf - AM * BM * pf;
   pb = pb - BM * AM * pb;

   %  Update forward and backward predictor coefficients - Eqs.(15.84),(15.85)
   if m ~= 1
      for k = 1:m-1
         temp1 = A(:,:,k);
         A(:,:,k)   = A(:,:,k)   + AM * B(:,:,m-k);
         B(:,:,m-k) = B(:,:,m-k) + BM * temp1;
      end
   end
   
   %  Update residuals
   Tef = ef;
   ef(1:NUMCHS,N:-1:m+1) = ef(:,N:-1:m+1) + AM * eb(:,N-1:-1:m);
   eb(1:NUMCHS,N:-1:m+1) = eb(:,N-1:-1:m) + BM * Tef(:,N:-1:m+1);

   %  Verify if model order is adequate
   if m == IP, A = -A; B = -B; break, end
end
