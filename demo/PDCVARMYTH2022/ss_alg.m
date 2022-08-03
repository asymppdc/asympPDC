%% SS_ALG
%      Calculate the spectral density matrix (SS) from A matrix and residues.
%
%% Syntax
%       SS = SS_ALG(A, e_cov, nFreqs)
%
%% Input arguments
%       A       - autoregressive coefficients matrix
%       e_cov   - residues
%       nFreqs  - number of frequencies
%
%% Output argument
%       SS      - Spectral density matrix
%
%   See also SS_ALG2, SS_ALG_B, SS_ALG_AB.
%

%   Copyright 2022 Koichi Sameshima and Luiz A. Baccala.


function SS = ss_alg(A, e_cov, nFreqs)

[nChannels, ~, ~] = size(A);
AL = A_to_f(A, nFreqs);
ss = zeros(size(AL));
for i = 1:nFreqs
   H = inv(reshape(AL(i,:,:),nChannels,nChannels));
   ss(i,:,:) = H * e_cov * H';
end

SS = permute(ss,[2,3,1]);
