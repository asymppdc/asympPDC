%% SS_ALG2
%     Calculate the two-sided frequency spectral density matrix (SS)
%
%% Syntax
%       SS = SS_ALG2(A, e_cov, nFreqs)
%
%% Input arguments
%       A       - autoregressive coefficients matrix
%       e_cov   - residues
%       nFreqs  - number of frequencies
%
%% Output argument
%       SS      - (nChannels, nChannels, 2*nFreqs) two-sided frequency
%                 spectral density matrix
%
% See also SS_ALG, SS_ALG_AB, SS_ALG_B
%

function SS = ss_alg2(A, e_cov, nFreqs)

[nChannels,~,~] = size(A);
AL = A_to_f2(A, nFreqs);

ss = zeros(size(AL));
for i = 1:2*nFreqs
   H = inv(reshape(AL(i,:,:),nChannels,nChannels));
   ss(i,:,:) = H * e_cov * H';
end;

SS = permute(ss,[2,3,1]);
