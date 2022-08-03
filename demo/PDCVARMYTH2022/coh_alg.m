%% COH_ALG
%       Calculate spectral coherence from power spectra, SS.
%
%% Syntax
%       Coh = COH_ALG(SS)
%
%% Input argument
%       SS      - Spectral density matrix
%
%% Output argument
%       Coh     - Squared spectral coherence
%

function Coh = coh_alg(SS)

[nChannels,nChannels2,nFreqs] = size(SS);
if nChannels == nChannels2, clear nChannels2; else error('Wrong SS dimension.'); end

Coh = zeros(size(SS));

% for k = 1:nFreqs,
for iu = 1:nChannels
   for ju = 1:nChannels
      Coh(iu,ju,:) = SS(iu,ju,:)./sqrt(SS(iu,iu,:).*SS(ju,ju,:));
   end
end
% end
