%% PAIRWISE_SPEC2
%       Calculate the two-sided frequency spectral density matrix (SS)
%
%% Syntax
%       Sx=PAIRWISE_SPEC2(u,nFreqs)
%
%% Input arguments
%       u       - (nChannels, data points) data sample
%       nFreqs  - number of frequencies
%
%% Output argument
%       Sx      - (nChannels x nChannels x 2*nFreqs) two-sided frequency
%                 spectral density matrix
%
%        See also SS_ALG2.
%          | <ss_alg2.html> |
%

function Sx=pairwise_spec2(u,nFreqs)

[nChannels,Ndata] = size(u);

SS = zeros(2,2,nFreqs); 
Sx = zeros(nChannels,nChannels,2*nFreqs);

for i = 1:nChannels
   y = u(i,:);
   [IP(i,i),pff,AA]  =  mvar(y,10,1,1);
   SS = ss_alg2(AA,pff,nFreqs);
   Sx(i,i,:,1) = SS;
end

k = 1;
for i = 1:nChannels
   for j = i+1:nChannels
      y = u([i j],:);
      [IP(i,j),pff,AA]  =  mvar(y,2*IP(i,i),1,5);
      k = k+1;
      SS = ss_alg2(AA,pff,nFreqs);
      Cx(i,j,:) = SS(1,2,:)./sqrt(SS(1,1,:).*SS(2,2,:));
      Sx(i,i,:,k) = SS(1,1,:);
      Sx(j,j,:,k) = SS(2,2,:);
   end
end

% Autospectra calculation
for i = 1:nChannels
   Su(i,i,:) = mean(Sx(i,i,:,:),4);
end

% Cross-spectra calculation
for i = 1:nChannels
   for j = i+1:nChannels
      Su(i,j,:) = Cx(i,j,:).*sqrt(Su(i,i,:).*Su(j,j,:));
      Su(j,i,:) = conj(Su(i,j,:));
   end
end
Sx = Su;
