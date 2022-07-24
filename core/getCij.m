function c = getCij(C,r,s,nFreq)
%function c = getCij(C,i,j,nFreq)
% 
% Input:      C [NumChannel, NumChannel, nFreqs], either PDC,Lpatnaik
%               Lv2inf, Lv2sup
%             [r,s] index-pair
%             nFreq - optional number of frequency values considered. 
%                     if not declared all freq values, i.e., nFreqs.
% Output:     c - C[r,s] element

if nargin < 4,
   [Nch Nch nFreq] = size(C);
end;
c=reshape(C(r,s,1:nFreq), nFreq,1,1,1,1);
