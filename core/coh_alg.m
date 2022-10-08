function Coh = coh_alg(SS)
%COH_ALG  Calculate spectral coherence from spectral density function.
%
% Syntax
%       Coh = coh_alg(SS)
%
% Input argument:
%        SS - (nChannels x nChannels x nFreqs) Power spectral density
%
% Output argument:
%       Coh - Spectral coherence

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.


   [m,n,nFreqs] = size(SS);
   Coh = zeros(size(SS));
   if m == n, nChannels = m; else error('Wrong SS dimension.'); end

   for k = 1:nFreqs
       for iu = 1:nChannels
           for ju = 1:nChannels
               Coh(iu,ju,k) = SS(iu,ju,k)./sqrt(SS(iu,iu,k).*SS(ju,ju,k));
           end
       end
   end
