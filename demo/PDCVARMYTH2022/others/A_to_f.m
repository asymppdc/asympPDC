%% A_to_f
%        Calculates A(f), in the positive frequency domain.
%
%% Syntax
%       AL = A_to_f(A, nFreqs)
%
%% Input arguments
%         A      - (nChannels x nChannels x p) Recurrence matrix 
%                  (nChannels - number of signals, p - model order)
%         nFreqs - frequency resolution
%
%% Output arguments
%         AL    - (nFreqs, nChannels, nChannels) A(f)
%
%        See also A_TO_F2.
%          | <A_to_f2.html> |
%

function AL = A_to_f(A, nFreqs)

[nChannels,~,p] = size(A);
Jimag = sqrt(-1);

% 'exponents' is an array of FFT exponents, on all frequency range for each lag.
exponents = reshape((-Jimag*pi*kron(0:(nFreqs-1),(1:p))/nFreqs),p,nFreqs).';

Areshaped = reshape(A, nChannels,nChannels,1,p);
Af = zeros(nChannels,nChannels,nFreqs,p);
for kk = 1:nFreqs
   Af(:,:,kk,:) = Areshaped;
end

for i = 1:nChannels
   for k = 1:nChannels
      Af(i,k,:,:) = reshape(Af(i,k,:,:),nFreqs,p).*exp(exponents);
   end
end

Af = permute(Af, [3,1,2,4]);
AL = zeros(nFreqs,nChannels,nChannels);
for kk = 1:nFreqs
   temp = zeros(nChannels,nChannels);
   for k = 1:p
      temp = temp+reshape(Af(kk,:,:,k),nChannels,nChannels);
   end
   temp = eye(nChannels)-temp;
   AL(kk,:,:) = reshape(temp,1,nChannels,nChannels);
end
