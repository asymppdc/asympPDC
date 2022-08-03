function [AA,BB,we,pe]=varmacomp(y,p,q);
% capsule
[pfx,A,~,~,efx,~] = mcarns(y,floor(50));
[AA,BB,we,pe] = varmalset(y,efx,p,q);