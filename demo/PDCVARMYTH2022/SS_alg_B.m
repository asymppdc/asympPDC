%% SS_alg_B
%       Calculate the spectral density matrix (SS), B(f) and spectral coherence
%       from VMA representation B matrix.
%
%% Syntax
%       [SS,VT,Coh]=SS_alg_B(B,pf,nFreqs,Ndata,flgNoCoh)
%
%% Input arguments
%        B      - (nChannels x nChannels x q+1) VMA repesentation 
%        pf     - white input covariance matrix  
%        nFreqs - number of desired frequency points
%        Ndata  - data length
%        flgNoCoh - any value if complex coherence is not desired.
%
%% Output arguments
%        SS     - spectral density matrix
%        VT     - B(f) --- frequency domain representation of B
%        Coh    - complex coherence
%
%% Used to compute: 
%      DTF:
%        c = wasymp_dtf(u,VT,pf,nFreqs,'diag',0,SS);
%  and
%      PDC: 
%        c = wasymp_pdc(u,VT,pf,nFreqs,'diag',0,SS);
%   
% See also SS_ALG, SS_ALG2, SS_ALG_AB
%          | <ss_alg.html> | <ss_alg2.html> | <ss_alg_AB.html> |
%

function [SS,VT,Coh]=SS_alg_B(B,pf,nFreqs,Ndata,flgNoCoh)

[nChannels,~,qt] = size(B);
mt = nChannels * nChannels;

BB = reshape(permute(B,[3 1 2]),qt,mt);
BB = [BB; zeros(2*nFreqs-qt,mt)];
fb = fft(BB); fb = fb(1:nFreqs,:);

if nargin < 5 % No coherence calculation
   Coh = zeros(nChannels,nChannels,nFreqs);
else
   Coh = [];
end

SS = zeros(nChannels,nChannels,nFreqs);

for i = 1:nFreqs
   V = reshape(fb(i,:),nChannels,nChannels);
   VT(:,:,i) = V;
   SS(:,:,i) = V * pf * V';
   ST = SS(:,:,i);
   if nargin == 5
      U = real(diag(ST)); U = 1./sqrt(U); U = diag(U);
      Coh(:,:,i) = U * ST * U;
   end
end

if Ndata~=1
   SS = 2*pi*SS/Ndata;
end
