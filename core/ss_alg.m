%% SS_ALG
%        Calculate spectral coherence from spectral density function.
%
%% Syntax:
%       SS = SS_ALG(A, e_cov, nf, nfcalc)
%
%% Input arguments:
%        A      - autoregressive matrix
%        e_cov  - residues
%        nf     - number of frequencies
%        nfcalc - fraction o frequency range to calculate 
%
%% Output argument:
%        SS     - (nChannels x nChannels x nfcalc) Power spectral density
%
%% See also: COH_ALG

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.


function SS = ss_alg(A, e_cov, nf, nfcalc)

if nargin < 4, nfcalc = nf; end
[~, nChannels,~] = size(A);
AL = A_to_f(A, nf);
ss = zeros(size(AL));
for i = 1:nfcalc
   H = inv(reshape(AL(i,:,:),nChannels,nChannels));
   ss(i,:,:) = H * e_cov * H';
end

SS=permute(ss,[2,3,1]);
end