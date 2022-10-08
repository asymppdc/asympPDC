function c = getCij(C,r,s,nFreq)
%GETCIJ  Retrive C[r,s] element from 3D martrix C [NumChannel, NumChannel, nFreqs]
%
% Syntax:
%        c = getCij(C,r,s,nFreq)
%
% Input arguments:
%               C [NumChannel, NumChannel, nFreqs], either PDC,Lpatnaik
%               Lv2inf, Lv2sup
%             [r,s] index-pair
%             nFreq - optional number of frequency values considered. 
%                     if not declared all freq values, i.e., nFreqs.
%
% Output argument:
%      c - C[r,s] element

% (C) Koichi Sameshima & Luiz A. Baccala, 2022. 
% See file license.txt in installation directory for licensing terms.


   if nargin < 4
      [~,~,nFreq] = size(C);
   end

   c=reshape(C(r,s,1:nFreq), nFreq,1,1,1,1);
   end