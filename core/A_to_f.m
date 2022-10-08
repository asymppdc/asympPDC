function AL = A_to_f(A, nf)
%A_to_f  Calculates A(f), A in the frequency domain.
%
% Syntax:
%       AL = A_to_f(A, nf)
%
% Input arguments:
%         A     - (nChannels x nChannels x p) Recurrence matrix,
%                  where nChannels is the number of signals, and p the model
%                  order.
%         nf    - frequency resolution
%
% Output argument:
%         AL    - (nf, nChannels, nChannels) A(f)

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.


   [nChannels, ~, p] = size(A); % Dummy variable ~ is equivalent to nChannels
   Jimag = sqrt(-1);

   % Variable 'exponents' is an array of FFT exponents, on all frequency
   % range for each lag.
   exponents = reshape((-Jimag*pi*kron(0:(nf-1),(1:p))/nf),p,nf).';

   % Af multiplies the exp(ar) by the matrix A, for all frequencies, the reshape
   % and transpose functions are tricks to make the vector calculation possible.

   Areshaped = reshape(A, nChannels,nChannels,1,p);

   Af = zeros(nChannels,nChannels,nf,p);
   for kk = 1:nf
      Af(:,:,kk,:) = Areshaped;
   end

   for i = 1:nChannels
      for k = 1:nChannels
         Af(i,k,:,:) = reshape(Af(i,k,:,:),nf,p).*exp(exponents);
      end
   end

   Af = permute(Af, [3,1,2,4]);

   AL=zeros(nf,nChannels,nChannels);

   for kk = 1:nf
      temp = zeros(nChannels,nChannels);
      for k = 1:p
         temp = temp+reshape(Af(kk,:,:,k),nChannels,nChannels);
      end
      temp = eye(nChannels)-temp;
      AL(kk,:,:) = reshape(temp,1,nChannels,nChannels);
   end
