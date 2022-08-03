%% PDC_TOT_P
%        Calculate total information PDC as proposed in [1].
%
%% Syntax
%        [pdct,pdc,pdcr,pdcp,spdc,y0i]=PDC_TOT_P(cipdc,pf)
%
%% Input arguments
%        cipdc  - Complex ipdc (nChannels x nChannels x nFreqs) 
%        pf     - Innovations covariance matrix (nChannels x nChannels)
%
%% Output arguments
%        pdct - total iPDC (complex)
%        pdc  - |iPDC|^2  (real)
%        pdcr - instantaneous terms (to remove in the future?)
%        pdcp - iiPDC (complex)
%        spdc - complex ipdc (nChannels x nChannels x nFreqs)
%        y0i  - second component of Residual directed iPDC in Eq. (27) [1]
%
%% References:
%    [1] Baccala LA, Sameshima K (2021). Frequency domain repercussions of
%        instantaneous Granger causality. Entropy 23(8), 1037
%                      <https://doi.org/10.3390/e23081037>
%
%            (This link may not work from within MATLAB Web browser. The
%            work-around is to copy the link into your favorite browser.)
%
%        See also ASYMP_PDC.
%               | <asymp_pdc.html> |
%

% This code segment will be incorporate into the 'asymp_pdc.m' function of
% Asymp_PDC package.

function [pdct,pdc,pdcr,pdcp,spdc,y0i]=pdc_tot_p(cipdc,pf)

[nChannels,~,nFreqs] = size(cipdc);
[ipf,g0] = dediag(pf,1);
g0 = diag(g0);
rho_I = ipf-eye(nChannels);
y0 = cipdc;
for ff = 1:nFreqs
    y01(:,:,ff) = g0 * y0(:,:,ff);
    yi(:,:,ff)  = rho_I * y01(:,:,ff);
    ye(:,:,ff)  = ipf * y01(:,:,ff);
end
spdc = y01;
y0i  = yi;
yi   = yi.*conj(y01);
y0   = conj(y01).*y01;                                           
y    = y0 + yi;
pdct = y;
pdc  = y0;
pdcr = yi;
pdcp = ye;
spdc = cipdc;
