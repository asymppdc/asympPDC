%% A_to_f2
%        Calculates A(f), in two-sided frequency domain.
%
%% Syntax
%       AL = A_to_f2(A, nFreqs)
%
%% Input arguments
%       A     - (nChannels x nChannels x p) Recurrence matrix 
%                  (nChannels - number of signals, p - model order)
%       nFreqs    - frequency resolution
%
%% Output arguments
%       AL    - (2*nFreqs, nChannels, nChannels) two-sided A(f)
%
%        See also A_TO_F.
%

function AL = A_to_f2(A, nFreqs)

[nChannels,~,p] = size(A);
Jimag = sqrt(-1);

% 'exponents' is an array of FFT exponents, on all frequency range for each lag.
exponents = reshape((-Jimag*pi*kron(0:(2*nFreqs)-1,(1:p))/nFreqs),p,2*nFreqs).';

Areshaped = reshape(A, nChannels,nChannels,1,p);
Af = zeros(nChannels,nChannels,2*nFreqs,p);
for kk = 1:2*nFreqs
   Af(:,:,kk,:) = Areshaped;
end

for i = 1:nChannels
   for k = 1:nChannels
      Af(i,k,:,:) = reshape(Af(i,k,:,:),2*nFreqs,p).*exp(exponents);
   end
end

Af = permute(Af, [3,1,2,4]);
AL = zeros(2*nFreqs,nChannels,nChannels);
for kk = 1:2*nFreqs
   temp = zeros(nChannels,nChannels);
   for k = 1:p
      temp = temp+reshape(Af(kk,:,:,k),nChannels,nChannels);
   end
   temp = eye(nChannels)-temp;
   AL(kk,:,:) = reshape(temp,1,nChannels,nChannels);
end
