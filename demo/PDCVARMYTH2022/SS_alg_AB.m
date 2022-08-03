%% SS_alg_AB
%       Calculate the spectral density matrix (SS), B(f) and spectral coherence
%       from VARMA representation A and B matrices.
%
%% Syntax 
%       [SS,VT,Coh]=SS_alg_AB(A,B,pf,nFreqs,Ndata,flgNoCoh)
%
%% Input arguments
%       A         - [nChannels,nChannels,p+1] - VAR repesentation part
%       B         - [nChannels,nChannels,q+1] - VMA repesentation part
%       pf        - white input covariance matrix  
%       nFreqs    - number of desired frequency points
%       Ndata     - data length
%       flgNoCoh  - any value if complex coherence is not desired.
%
%% Output arguments
%       SS        - spectral density matrix
%       VT        - B(f) frequency domain representation of B
%       Coh       - complex coherence
%
%% Used 
%       To compute dtf:
%          c=wasymp_dtf(u,VT,pf,nFreqs,'diag',0,SS);
%       To compute pdc: 
%          c=wasymp_pdc(u,VT,pf,nFreqs,'diag',0,SS);
%
% See also SS_ALG, SS_ALG2, SS_ALG_B
%          | <ss_alg.html> | <ss_alg2.html> | <ss_alg_B.html> |
%

function [SS,VT,Coh]=SS_alg_AB(A,B,pf,nFreqs,Ndata,flgNoCoh)

[nChannels,~,pt]  = size(A);
[mt, ~, qt] = size(B);

if nChannels ~= mt, error('A and B have different numbers of channel.'); end
mt = nChannels * nChannels;

AA = reshape(permute(A,[3 1 2]),pt,mt);
AA = [reshape(eye(nChannels),1,mt); -AA; zeros(2*nFreqs-pt-1,mt)];
fa = fft(AA); fa = fa(1:nFreqs,:);

BB = reshape(permute(B,[3 1 2]),qt,mt);
BB = [BB; zeros(2*nFreqs-qt,mt)];
fb = fft(BB); fb = fb(1:nFreqs,:);

if nargin < 6
   Coh = zeros(nChannels,nChannels,nFreqs);
else
   Coh = [];
end

SS = zeros(nChannels,nChannels,nFreqs);

for i = 1:nFreqs
   V  = reshape(fb(i,:),nChannels,nChannels);
   VI = reshape(fa(i,:),nChannels,nChannels);
   V  = inv(VI) * V;
   VT(:,:,i) = V;
   SS(:,:,i) = V * pf * V';
   ST = SS(:,:,i);
   if nargin < 6
      U = real(diag(ST));
      U = 1./sqrt(U);
      U = diag(U);
      Coh(:,:,i) = U * ST * U;
   end
end

if Ndata ~= 1
   SS = 2*pi*SS/Ndata;
end
