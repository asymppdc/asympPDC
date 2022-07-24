%% SS_ALG
%        Calculate spectral coherence from spectral density function.
%% Syntax
%       SS = SS_ALG(A, e_cov, nf, nfcalc)
%% Input arguments
%        A      - autoregressive matrix
%        e_cov  - residues
%        nf     - number of frequencies
%        nfcalc - fraction o frequency range to calculate 
%% Output arguments
%        SS     - (nChannels x nChannels x nfcalc) Power spectral density
%

function SS = ss_alg(A, e_cov, nf, nfcalc)

if nargin < 4, nfcalc = nf; end;
[nnChannels, nChannels, r] = size(A);
AL = A_to_f(A, nf);
ss = zeros(size(AL));
for i = 1:nfcalc,
   H = inv(reshape(AL(i,:,:),nChannels,nChannels));
   ss(i,:,:) = H * e_cov * H';
end;

SS=permute(ss,[2,3,1]);
