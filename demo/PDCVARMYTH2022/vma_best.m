%% VMA_BEST
%       compute vector MA model
%% Syntax
%       [IP,pf,B,vaic,Vaicv] = VMA_BEST(u,maxIP,criterion)
%% Input arguments
%       u         - signal Nchannels x uration
%       MaxIP     - maximum model order (>0)
%       criterion - model order choice
%
%% Output argumentS
%       IP      - (q) for MA(q) value
%       pf      - Innovations covariance matrix
%       B       - Matrix from model,
%       vaic    - AIC value
%       Vaicv   - Vector evolution of aic values
%

function [IP,pf,B,vaic,Vaicv] = vma_best(u,maxIP,criterion)

[nChannels,ESL] = size(u);
% [pfx,A,~,~,efx,~] = mcarns(u,floor(.4*length(u)));
[pfx,A,~,~,efx,~] = mcarns(u,50);

i = 0;
[va,ve,vpf] = vmalse(u,efx,i);
npf = vpf;

% criterion routines changed to suit number of parameters
switch criterion
   case 1,  % Akaike's Informaion Criterin (AIC)
      vaic = ESL*log(det(npf)) ...
                         + 2*(nChannels*nChannels)*(i+1); %(4.3.2)(Lutkepohl'05)
   case 2,  % Hannan-Quinn (HQ)
      vaic = ESL*log(det(npf)) + 2*log(log(ESL)) ...
                          *(nChannels*nChannels)*(i+1);   %(4.3.8)(Lutkepohl'05)
   case 3,  % Schwarz (SC) (Schwarz, 1978)
      vaic = ESL*log(det(npf)) ...
                  + log(ESL)*(nChannels*nChannels)*(i+1); %(4.3.9)(Lutkepohl'05)
   case 4,  % FPE - Final prediction error (Akaike, 1970)
      vaic = log(det(npf)*((ESL+nChannels*(i+1)+1) ...
                    /(ESL-nChannels*(i+1)-1))^nChannels); %(4.3.1)(Lutkepohl'05)
end

Vaicv(1) = vaic;
naic = vaic;

for i = 1:maxIP
   [na,ne,npf] = vmalse(u,efx,i);
   switch criterion
      case 1,  % Akaike's Informaion Criterin (AIC)
         naic = ESL*log(det(npf)) ...
                    + 2*(nChannels*nChannels)*(i+1);      %(4.3.2)(Lutkepohl'05)
      case 2,  % Hannan-Quinn (HQ)
         naic = ESL*log(det(npf)) + 2*log(log(ESL)) ...
                            *(nChannels*nChannels)*(i+1); %(4.3.8)(Lutkepohl'05)
      case 3,  % Schwarz (SC) (Schwarz, 1978)
         naic = ESL*log(det(npf)) + log(ESL) ...
                            *(nChannels*nChannels)*(i+1); %(4.3.9)(Lutkepohl'05)
      case 4,  % FPE - Final prediction error (Akaike, 1970)
         naic = log(det(npf)*((ESL+nChannels*(i+1)+1) ...
                 /(ESL-nChannels*(i+1)-1))^nChannels);    %(4.3.1)(Lutkepohl'05)
         %       otherwise
         %nop
   end

   Vaicv(i+1) = naic;
   if vaic > naic
      vaic = naic;
      va = na;
      vpf = npf;
      ve = ne;
   else
      IP = i-1;
      B = va;
      ef = ve;
      pf = pfx;
      return
   end
   IP = i-1;
   B = va;
   ef = ve;
   pf = pfx;
end