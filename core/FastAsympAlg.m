%% FASTASYMPALG
%       Modified Fast asymptotic PDC and DTF calculation algorithm implemented
%       by Rezaei et al. (2022) to which c struct variable is added with fields
%       that make it compatible with xplot PDC and DTF pretty plotting routine
%       of the asympPDC Package. Additionally some syntax compatibility issues
%       to older MATLAB versions, previous to R2020b, were solved as the
%       original FastAsympAlg function seems to be implemented using the most
%       recent version of MATLAB (after R2020a).
%
%% Syntax
%       [c,d] = FASTASYMPALG(s,A,pf,Freq,measure,SR,alpha)
%
%% Input Arguments:
%       s      : signal - [time x nchannel]
%       A      : Matrix of coefficients estimated by MVAR
%       pf     : Covariance of noise
%       SR	   : Sample Rate
%       alpha  : significance level of statistical test
%       Freq   :  frequency band
%       measure: 'DTF, iDTF, gDTF|DC ( with metric='gamma_group' )
%                'PDC, iPDC, gPDC (with metric 'pi_group')
%
%% Output Arguments (added for compatibility with asympPDC Package)
%        
%        c: Original struct variable with the following fields:
%        |-- .Phi       - |PDC|^2 or |DTF|^2 estimates
%        |-- .Pval      - p-values associated to PDC2 or DTF2 estimates. 
%        |-- .Threshold - |PDC|^2 threshold value with (1-avalue) significance 
%        |                level.
%        |-- .Phi_th    - (alpha-1)% statistically significant Connectivity 
%        |                measures
%        |-- .CIupperbound - upper and lower (1 - alpha) confidence interval of 
%        |-- .CIlowerbound   |PDC|^2 or |DTF|^2 estimates
%        |-- .measure   - Type of Phi: 'pdc','gpdc','ipdc','dtf','dc' or 'idtf' 
%        |-- .SS        - Spectral density
%        |-- .Comment   - Comment
%      i.e.
%        .{Phi,Pval,Threshold,Phi_th,CIupperbound,CIlowerbound,measure,SS,
%                                                                       Comment}
%
%        d: asympPDC Package compatible struct variable containing following
%        |  fields used by xplot and xplot_pvalues routines:
%        |-- .pdc2/dtf2 - |PDC|^2 / |DTF|^2  estimates
%        |-- .pvalues   - p-values associated to PDC2 / DTF2  estimates. 
%        |-- .th        - |PDC|^2 / |DTF|^2  threshold value with (1-avalue) 
%        |                significance level.
%        |-- .pdc2_th/dtf2_th - significant values above threshold level
%        |                otherwise NaN in freqs
%        |-- .{ci1,ci2} - upper and lower (1 - alpha) confidence interval of 
%        |                |PDC|^2 / |DTF|^2  estimates
%        |-- .metric    - metric used for |PDC|^2 or |DTF|^2  estimation:
%        |                'euc', 'diag' or 'info'
%        |-- .alpha     - significance level
%        |-- .p         - VAR model order
%        |-- .Comment   - comment
%        |-- .SS        - power spectra
%        +-- .coh2      - squared spectral coherence
%      i.e.
%        .{{pdc2|dtf2},pvalues,th,{pdc2_th|dtf2_th},ci1,ci2,metric,alpha,p,
%                                                               Comment,SS,coh2}
%% See also: FASTASYMPALG_ORIG, ASYMP_PDC, ASYMP_DTF
%

%% References (added for completeness of references)
%
% [1] Kaminski, M.J. and Blinowska, K.J., 1991. A new method of the description 
%     of the information flow in the brain structures. Biol Cybern 65, 
%     pp.203-210,
%     <https://doi.org/10.1007/bf00198091>
%
% [2] Baccalá, L.A and Sameshima, K., 2001. Partial directed coherence: A new
%     concept in neural structure determination. Biol Cybern 84, pp.463-474.
%     <https://doi.org/10.1007/PL00007990>
%
% [3] Takahashi, D.Y., Baccalá, L.A., and Sameshima, K., 2007. Connectivity
%     inference between neural structures via partial directed coherence. J
%     Appl Stat 34(10), pp.1259-1273
%     <https://doi.org/10.1080/02664760701593065>
%
% [4] Baccalá, L.A., De Brito, C.S., Takahashi, D.Y. and Sameshima, K., 2013.
%     Unified asymptotic theory for all partial directed coherence forms. Phil
%     Trans R Soc A, 371(1997), p.20120158.
%     <https://doi.org/10.1109/TBME.2016.2550199>
%
% [5] Baccalá, L.A., Takahashi, D.Y. and Sameshima, K., 2016. Directed transfer
%     function: Unified asymptotic theory and some of its implications. IEEE
%     Trans Biomed Eng, 63(12), pp.2450-2460.
%     <https://doi.org/10.1109/TBME.2016.2550199>
%
% [6] Rezaei, F., Alamoudi, O., Davani, S. and Hou, S., 2022. Fast
%     asymptotic algorithm for real-time causal connectivity analysis of
%     multivariate systems and signals. Signal Process, 204, 108822.
%     <https://doi.org/10.1016/j.sigpro.2022.108822>
%
%% Written by: 
%    Farnaz Rezaei  2019/04/10 
%    - Last version 2020/03/07
%    - Obtained from <https://github.com/Farnaz-codes/FAA> 2022/11/22
%
%    Adapted and modified to turn it compatible with previous MATLAB releases,
%    compatible down to MATLAB R2014b and Octave versions 6.2, 6.3, 7.1, 7.2 and
%    7.3 as partially tested on 2022/12/10 
%------------------------------------------------------------------------------%

function [c,d] = FastAsympAlg(s,ARcoef,Ecov,Freq,measure,SR,alpha)

% ------------------------------------------------------------------------------
% This function was modified/adapted to deal with rdivide compatibility issue
% for older version of MATLAB
% ------------------------------------------------------------------------------

if ~isOctave()
   v = ver('MATLAB');
   versionNumberMatlab = str2double(v.Version);
else
   versionNumberMatlab = 8.4; % Does Octave have same compatibility as ML R2014b (?!)
end
% For MATLAB release history see at <https://en.wikipedia.org/wiki/MATLAB>
% Release      Version
%  2014b         8.4
%  2017a         9.2
%  2020a         9.9
%  2021b         9.11

% ------------------------------------------------------------------------------
% Checking input arguments
% ------------------------------------------------------------------------------

if (nargin<7 || alpha==0)
   alpha = 0;
   c.Comment = ['Estimation ', upper(measure)];
else
   c.Comment = ['Asymptotic distribution of ', upper(measure), ' and its ', ...
                               num2str(100*(1-alpha)),' % Confidence interval'];
end

if isempty(SR); SR = 2*Freq(end); end

measure = lower(measure);
c.measure = measure;
metric = 'pi_group';
if isempty(regexp(measure,'pdc', 'once')); metric = 'gamma_group'; end
if (isempty(regexp(measure,'pdc', 'once')) && ...
                  isempty(regexp(measure,'dtf', 'once')) || length(measure) > 4)
   error('Unknown metric.')
end
if (length(measure)==4 && (strcmp(measure(1),'i')==0 && ...
                                                     strcmp(measure(1),'g')==0))
   error('Unknown metric.')
end
[np, nChannels] = size(s);
if np < nChannels, s = s'; [np,nChannels] = size(s); end

% ------------------------------------------------------------------------------
% Noncompatible syntax for MATLAB older than R2016b for
% s = s-mean(s); 
%
if versionNumberMatlab < 9.2  % MATLAB previous to R2017a
   s = bsxfun(@minus,s, mean(s)); % 2022.11.07  KS
   % ---------------------------------------------------------------------------
   % Note: Workaround without using bsxfun() could be [2022.11.07 KS]
   % ===============================
   % [~, kanais]=size(s);
   % smean = mean(s);
   % for k = 1:kanais,
   %    s(:,k) = s(:,k) - smean(k);
   % end
   % ===============================
else
   s = s - mean(s);
end

icdf_norm_alpha = norminv(1-alpha/2.0,0,1);
if ~issymmetric(Ecov)
   Ecov = triu(Ecov,1) + triu(Ecov)';
end
%------------------------------------------------------------------------------%
p = size(ARcoef,3);                   % model order
AR = reshape(ARcoef,nChannels,nChannels*p);
Ecovdotproduct = Ecov.*Ecov';     % Hadamard product of Covariance of noise

% ------------------------------------------------------------------------------
%                            Sn and Sd table(4)
%
switch metric
   case 'gamma_group'
      evar_d = eye(nChannels).*Ecov;
      if strcmp(measure,'dtf');evar_d = eye(nChannels);end
      Sd = evar_d;
      if strcmp(measure,'idtf');Sd = Ecov;end
      
   case 'pi_group'
      evar_d = inv((eye(nChannels).*Ecov));
      if strcmp(measure,'pdc');evar_d = eye(nChannels);end
      Sd = evar_d;
      if strcmp(measure,'ipdc');Sd = pinv(Ecov);end
      ce = diag(evar_d).*diag(Ecov);
end

% ------------------------------------------------------------------------------
%                       Variables initialization
%
Phi     = zeros(nChannels, nChannels, length(Freq));
pvalues = zeros(nChannels, nChannels, length(Freq));
TH      = zeros(nChannels, nChannels, length(Freq));
CIup    = zeros(nChannels, nChannels, length(Freq));
CIlow   = zeros(nChannels, nChannels, length(Freq));
SS      = zeros(nChannels, nChannels, length(Freq));
gamma   = zeros(nChannels*p, nChannels*p);

% ------------------------------------------------------------------------------
%       Estimation of gamma and its decomposition
%
for m = 1:p
   for n = m:p
      SlagOut = s(n-m+1:(end-m+1),:)'*s(1:(end-n+1),:)/np;
      gamma(((m-1)*nChannels+1):m*nChannels, ...
                                    ((n-1)*nChannels+1):n*nChannels) = SlagOut;
      gamma(((n-1)*nChannels+1):n*nChannels, ...
                                    ((m-1)*nChannels+1):m*nChannels) = SlagOut';
   end
end
invgamma = pinv(gamma);

% ------------------------------------------------------------------------------
%                          Frequency loop
%
for ff = 1:length(Freq)
   f1 = Freq(ff);
   f  = (f1-1)/(SR);  %f starting at 0
   C1 = cos(2*pi*f*(1:p));
   S1 = sin(2*pi*f*(1:p));
   C  = [C1;-S1];
   Af = eye(nChannels)-AR*kron(C1-1i*S1,eye(nChannels)).';
   %%-------------------------------------------------------------------------%
   Hf = pinv(Af);
   SS(:,:,ff) = Hf*Ecov*Hf';
   OmegaX    = zeros(nChannels,nChannels);
   OmegaEcov = zeros(nChannels,nChannels);
   % ---------------------------------------------------------------------------
   % Estimation of measures
   % phi = x'*Tn*x/x'*Td*x;
   switch metric
      case 'gamma_group'
         X = Hf;
         num = abs(X).^2*evar_d;
         den = sum(num,2);
         if strcmp(measure,'idtf')
            den = diag(real(X*Ecov*X'));
         end
      case 'pi_group'
         X = Af;
         num = evar_d*abs(X).^2;
         den = sum(num);
         if strcmp(measure,'ipdc')
            den = diag(real(X'*pinv(Ecov)*X));
            den = den';
         end
   end
   
   % ---------------------------------------------------------------------------
   % Note: A workaround to perform     2022.11.07 KS
   % %   Phi(:,:,ff) = num./den;
   %
   %     [~, kanais]=size(s);
   %     tmp = zeros(kanais,kanais);
   %     for k = 1:kanais,
   %        tmp(:,k) = num(:,k)/den(k)
   %     end
   %     Phi(:,:,ff) = tmp; % (N,N) ./ (N, 1) not compatible with ML R2013 
   
   %     Phi(:,:,ff) = num./den; % (N,N) ./ (N, 1) not compatible with ML R2017a
   
   if versionNumberMatlab < 9.2  % MATLAB before release 2016b
      Phi(:,:,ff) = bsxfun(@rdivide, num, den); % 2022.11.07 KS      
   else
      Phi(:,:,ff) = num./den; % (N,N) ./ (N, 1) not compatible with ML < R2016b
   end
   
   %-------------------------------------------------------------------------%
   if alpha~=0
      % If the Null hypothesis is rejected, Threshold is: (I add this if
      % later) FR
      switch metric
         case 'gamma_group'
            %--------------------------------------------------------------%
            % Estimating df (degree of freedom and multiplier)
            % without eigenvalues estimation and cholesky factors
            % Defining gamma-group parameters
            F1 = 0.5*[1,1j;-1j,1];
            HEcov = X*Ecov;
            miu1 = sum(HEcov.*X,2);
            miu2 = sum(HEcov.*conj(X),2);
            temp = kron(eye(p),X.')*invgamma;
            Gp_gamma_1 = temp*kron(eye(p),X);
            Gp_gamma_2 = temp*kron(eye(p),conj(X));
            %---------------------------------------------------------%
            T_th_1 = nan(1,nChannels); T_th_2 = T_th_1; T_th_3 = T_th_1;
            P2 = nan(2,2,nChannels); P1 = P2;
            for j = 1:nChannels
               tempj_1 = Gp_gamma_1(j:nChannels:nChannels*p, ...
                                                       j:nChannels:nChannels*p);
               tempj_2 = Gp_gamma_2(j:nChannels:nChannels*p, ...
                                                       j:nChannels:nChannels*p);
               P2(:,:,j) = F1*C*tempj_2*C'*F1;
               P1(:,:,j) = F1*C*tempj_1*C'*F1.';
               
               T_th_1(j) = real(trace(P2(:,:,j))); %trace (P2)
% ------------------------------------------------------------------------------
               if versionNumberMatlab <= 9.2 % MATLAB before release 2018a (?)
                  % T_th_2(j) = real(sum(P2(:,:,j).*P2(:,:,j).','all'));% trace (P2*P2)                                      ^^^
                  % T_th_3(j) = real(sum(P1(:,:,j).*P1(:,:,j)','all')); % trace (P2*P2)
                  %                                             ^^^
                  %   "S = sum(A,'all') computes the sum of all elements of A. 
                  % This syntax was introduced for MATLAB® versions R2018b
                  %    and newer. (www.mathworks.com)"
                  % A workaround is to use 'sum(sum(sum( argument )))' as
                  % the argument A has three dimensions. % 2022.11.07 KS 

                  T_th_2(j) = real(sum(sum(sum(P2(:,:,j).*P2(:,:,j).'))));                
                  T_th_3(j) = real(sum(sum(sum(P1(:,:,j).*P1(:,:,j)'))));
               else
                  T_th_2(j) = real(sum(P2(:,:,j).*P2(:,:,j).','all')); % trace (P2*P2)
                  T_th_3(j) = real(sum(P1(:,:,j).*P1(:,:,j)','all'));  % trace (P2*P2)                  
               end
% ------------------------------------------------------------------------------
            end
            mult = abs(miu1).^2 * T_th_3 + abs(miu2).^2 * T_th_2;
            patdf = (2*(abs(miu2)*T_th_1).^2) ./ mult;
            patmul = (mult*evar_d) ./ ((den.*abs(miu2))*T_th_1);
            %-----------------------------------------------------------------%
         case 'pi_group'
            cgamma = nan(1,nChannels); nugamma = nan(1,nChannels);
            for j = 1:nChannels
               temp = C*invgamma(j:nChannels:nChannels*p, ...
                                                    j:nChannels:nChannels*p)*C';
               cgamma(j)  = sum(temp(:).^2)/(den(j)*trace(temp));
               nugamma(j) = trace(temp).^2/sum(temp(:).^2);
            end
            patmul = ce*cgamma;
            patdf  = repmat(nugamma,nChannels,1);
      end
% ------------------------------------------------------------------------------
%     The icdf function was introduced before R2006a but is incompatible(?) with
%     Octave, and, additionally, the corresponding specific function is faster
%     than the generic one. For instance, icdf('chi2', ...) can be substituted
%     by chi2cdf( ) and icdf('Normal', ...) by norminv( ). % 2022.11.07 KS
%
      %TH(:, :, ff) = real(patmul).*icdf('chi2',(1-alpha), real(patdf))./(np);
      %                             ============
      % should be replaced by       vvvvvvvv
      TH(:, :, ff) = real(patmul).* chi2inv((1-alpha), real(patdf)) ./ np;
      %                             ========
      pvalues(:,:,ff) = 1 - chi2cdf(Phi(:,:,ff)*np./real(patmul),real(patdf));

% ------------------------------------------------------------------------------
      % Confidence interval
      % Covariance of measures with respect to the X :OmegaX
      switch metric
         case 'gamma_group'
            % -----------------------------------------------------------------%
            %Edited - 8/25/2020 FR
            Xbar = [real(X),imag(X)];
            temp = Xbar*kron(F1*C,Sd);
            gsigma1 = sum(temp*Gp_gamma_1.*temp,2);
            gsigma2 = sum(temp*Gp_gamma_2.*conj(temp),2);
            T2 = real(miu1.*gsigma1+miu2.*gsigma2);
            for j = 1:nChannels
               temp1 = sum(Xbar(:,[j,j+nChannels])*P1(:,:,j).*Xbar(:, ...
                                                            [j,j+nChannels]),2);
               temp2 = sum(Xbar(:,[j,j+nChannels])*P2(:,:,j).*Xbar(:, ...
                                                            [j,j+nChannels]),2);
               T1 = real(miu1.*temp1 + miu2.*temp2);
               temp3 = Xbar(:,[j,j+nChannels]) * F1 * C;
               temp3_1 = sum((temp3*Gp_gamma_1(j:nChannels:nChannels*p, ...
                                                                   :)).*temp,2);
               temp3_2 = sum((temp3*Gp_gamma_2(j:nChannels:nChannels*p, ...
                                                             :)).*conj(temp),2);
               T3 = 2*real(miu1.*temp3_1+miu2.*temp3_2);
               OmegaX(:,j) = 8*(evar_d(j,j)^2*T1 + ...
                      Phi(:,j,ff).^2.*T2-evar_d(j,j)*Phi(:,j,ff).*T3)./(den.^2);
            end
            %--------------------------------------------------------%
            % updated 6/28/2020 FR
         case 'pi_group'
            for j = 1:nChannels
               Gp = 4*C*invgamma(j:nChannels:nChannels*p, ...
                                         j:nChannels:nChannels*p)*C'/(den(j)^2);
               Xbar = [real(X(:,j)),imag(X(:,j))];
               temp1 = Xbar*Gp*Xbar';
               switch measure
                  case 'pdc'
                     temp2 = Ecov.*temp1;
                     T1 = diag(Ecov).*diag(temp1);
                     T2 = sum(temp2(:));
                     T3 = sum(temp2,2);
                     
                  case 'gpdc'
                     temp2 = Ecov.*temp1;
                     T1 = diag(evar_d).*diag(temp1);

% ------------------------------------------------------------------------------
%      The expression "T3 = diag(evar_d).*temp2*diag(evar_d)" requires
%      parenthesis to guarantee proper order of operations for older
%      versions of MATLAB
%
%                    % T3 = diag(evar_d).*temp2*diag(evar_d);
%                    %  Arrays dimension:
%                    %  size(diag(evar_d)) : (nChannels x 1)
%                    %  size(temp2)        : (nChannels x nChannels)
%                    %  size(diag(evar_d)) : (nChannels x 1)
%
%               % T3 =   diag(evar_d) .*         temp2         * diag(evar_d);
%               % T3 = (nChannels x 1).*(nChannels x nChannels)*(nChannels x 1)
%               %       ============= ^^ =====================
%               %      So this pair has incompatible array sizes. 
%               %      bsxfun with vvvvvvvv resolves the issue, apparently.
%                    T3   = bsxfun(@rdivide,diag(evar_d),temp2)*diag(evar_d);
%                    or
                     T3 = diag(evar_d).*(temp2*diag(evar_d));
%                using parentheses      ^                  ^                     
% ------------------------------------------------------------------------------
                     T2 = sum(T3);
                     
                  case 'ipdc'
                     temp2 = Sd.*temp1;
                     T1 = diag(evar_d).*diag(temp1);
                     T2 = sum(temp2(:));
                     T3 = T1;
               end
               OmegaX(:,j) = T1+Phi(:,j,ff).^2.*T2-2*Phi(:,j,ff).*T3;
            end
            %---------------------------------------------------------%
      end
      % end of Confidence interval - Part 1

% ------------------------------------------------------------------------------
      % 2) Covariance of measures with respect to the noise :OmegaEcov  
      %    updated 6/19/2020 FR
      if (~strcmp(measure,'pdc') || ~strcmp(measure,'dtf') )
         switch metric
            case 'pi_group'
               T1 = (evar_d^2)*(abs(X).^4);
               switch measure
                  case 'gpdc'
                     AbsB2 = abs(X).^2;
                     temp2 = AbsB2'*(evar_d^2)*Ecovdotproduct*(evar_d^2)*AbsB2;
                     T2 = repmat(diag(temp2)',nChannels,1).*Phi(:,:,ff).^2;
                     T3 = AbsB2.*Phi(:,:,ff).*((evar_d^2)*Ecovdotproduct* ...
                        (evar_d^2)*AbsB2);
                     
                     if versionNumberMatlab <= 9.2  % MATLAB release 2017a
                        OmegaEcov = bsxfun(@rdivide, 2*(T1+T2-2*T3), den.^2);  
                        
                     elseif versionNumberMatlab < 9.9 % MATLAB release 2020a
                        OmegaEcov = (2*(T1+T2-2*T3))./(den.^2);
                        %                             ^      ^  KS 2022.11.07
                     else
                        OmegaEcov = 2*(T1+T2-2*T3)./den.^2; % as implemented by
                                                            % Rezaei et al.(2022)
                     end

                     %--------------------------------------------------%
                  case 'ipdc'
                     T3 = T1.*Phi(:,:,ff);
                     temp2 = 0.5*(abs(X'*Sd*X).^2+abs(X.'*Sd*X).^2);
                     T2 = repmat(diag(temp2)',nChannels,1).*Phi(:,:,ff).^2;
                     
                     if versionNumberMatlab < 9.2  % MATLAB before release 2017a
                        OmegaEcov = bsxfun(@rdivide, 2*(T1+T2-2*T3), den.^2);  % 2022.11.07 KS
                        %           ================
                     elseif versionNumberMatlab < 9.9 % MATLAB before release 2020a
                        OmegaEcov = (2*(T1+T2-2*T3))./(den.^2);
                        %                             ^      ^  KS                        
                     else
                        OmegaEcov = 2*(T1-2*T3+T2)./den.^2; % as implemented by Farnaz Rezaei
                     end

               end
               %---------------------------------------------------------------%
            case 'gamma_group'
               T1 = abs(X).^4*(evar_d^2);
               AbsH2 = abs(X).^2;
               switch measure
                  case 'gdtf'
                     temp2 = diag(AbsH2*Ecovdotproduct*AbsH2');
                     T2 = repmat(temp2,1,nChannels).*Phi(:,:,ff).^2;
                     T3 = AbsH2.*Phi(:,:,ff).*(AbsH2*Ecovdotproduct);
                     
                     if versionNumberMatlab < 9.2  % MATLAB before release 2017a
                        OmegaEcov = bsxfun(@rdivide, 2*(T1+T2-2*T3), ...
                                   repmat(den,1,nChannels).^2);  % 2022.11.07 KS
                        
                     elseif versionNumberMatlab < 9.9 % MATLAB before release 2020a
                        OmegaEcov = (2*(T1+T2-2*T3))./(repmat(den,1,nChannels).^2);
                        % Using parentheses           ^                          ^
                     else
                        OmegaEcov = 2*(T1+T2-2*T3) ./ repmat(den,1,nChannels).^2;
                     %              ============== ^^ =========================
                     %             these are incompatible array sizes. 
                     %             bsxfun resolves the issue.                        
                     end
                     %--------------------------------------------------%
                  case 'idtf'
                     T3 = Phi(:,:,ff).*AbsH2.*(HEcov.*(conj(X)*Ecov));
                     temp2 = 0.5*(abs(miu1).^2+abs(miu2).^2);
                     T2 = repmat(temp2,1,nChannels).*(Phi(:,:,ff).^2);

                     if versionNumberMatlab < 9.2  % MATLAB before release 2017a
                        OmegaEcov = bsxfun(@rdivide, 2*(T1+T2-2*T3), ...
                                   repmat(den,1,nChannels).^2);  % 2022.11.07 KS
                        
                     elseif versionNumberMatlab < 9.9 % MATLAB before release 2020a
                        OmegaEcov = (2*(T1+T2-2*T3))./(repmat(den,1,nChannels).^2);
                        % Using parentheses           ^                          ^                        
                     else
                        OmegaEcov = 2*(T1+T2-2*T3) ./ repmat(den,1,nChannels).^2;
                        %           =============  ^^ ==========================
                        %           These are incompatible array sizes. 
                        %           bsxfun resolves the issue, apparently.                        
                     end

               end
               %---------------------------------------------------------------%
         end
      end   % end of Confidence interval-part2
      %------------------------------------------------------------------------%
      OmegaT = (OmegaEcov + OmegaX)/np;
      CIup(:,:,ff)  = Phi(:,:,ff) + sqrt(OmegaT)*icdf_norm_alpha;
      CIlow(:,:,ff) = Phi(:,:,ff) - sqrt(OmegaT)*icdf_norm_alpha;
   end
end
Phi_th = ((abs(Phi)-abs(TH)) > 0).*Phi;
Phi_th(Phi_th==0) = nan;
%-------------------------------------------------------------------------%
if alpha~=0
   c.CIupperbound = CIup;
   c.CIlowerbound = CIlow;
   c.Pval = pvalues;
   c.Phi = Phi;
   c.Comment = ['Asymptotic distribution of ', upper(measure), ' and its ', ...
                               num2str(100*(1-alpha)),' % Confidence interval'];
   c.Threshold = TH;
   c.Phi_th = Phi_th;
   c.SS = SS;
else
   c.Phi = Phi;
   c.SS = SS;
end

% ------------------------------------------------------------------------------
% Additional output struct variable compatible with xplot.m routine 
% in the asympPDC Package was added. KS
% (Single precision used to save memory space occupation)
%
if isempty(regexp(measure,'pdc', 'once'))
   metric = 'gamma_group'; 
   if ~isempty(regexp(measure,'idtf', 'once'))
      metric2 = 'info';
   elseif ~isempty(regexp(measure,'dc', 'once')) || ...
                                        ~isempty(regexp(measure,'gdtf', 'once'))
      metric2 = 'diag';
   else
      metric2 = 'euc';
   end
else
   if ~isempty(regexp(measure,'ipdc', 'once'))
      metric2 = 'info';
   elseif ~isempty(regexp(measure,'gpdc', 'once'))
      metric2 = 'diag';
   else
      metric2 = 'euc';
   end
end


if alpha ~= 0
   d.ci2 = single(CIup);
   d.ci1 = single(CIlow);
   d.pvalues = single(pvalues);
   d.metric = metric2;
   d.alpha = alpha;
   d.p = p;
   if strcmp(metric,'pi_group')
      d.pdc2 = single(Phi);
      d.pdc2_th = single(Phi_th);  % Significant PDC2
   else % strcmp(metric,'gamma_group')
      d.dtf2 = single(Phi);
      d.dtf2_th = single(Phi_th);  % Significant DTF2
   end
   d.Comment = ['Asymptotic distribution of ', upper(measure), ' and its ', ...
                sprintf('%5.3g',(100*(1-alpha))),' % Confidence interval'];
   d.th = single(TH);
   d.SS = single(SS);
   d.coh2 = single(coh_alg(SS)); % coh_alg function from KS+LAB asympPDC Package
else
   if strcmp(metric,'pi_group')
      d.pdc2 = single(Phi);
   else % strcmp(metric,'gamma_group')
      d.dtf2 = single(Phi);
   end
   d.alpha = alpha;
   d.p = p;
   d.SS = single(SS);
   d.metric = metric2;
   d.coh2 = single(coh_alg(SS)); % coh_alg function from KS+LAB asympPDC Package
end


%% Notes

%% TO DO LIST: 
% 1. add {pdc + dtf (Complex)} to d struct variable.

% ------------------------------------------------------------------------------
%% Output Arguments 
%        
%        c: original structure variable with the following fields:
%        |-- .Phi       - |PDC|^2 or |DTF|^2 estimates
%        |-- .Pval      - p-values associated to PDC2 or DTF2 estimates. 
%        |-- .Threshold - |PDC|^2 threshold value with (1-avalue) significance 
%        |                level.
%        |-- .Phi_th    - (alpha-1)% statistically significant Connectivity 
%        |                measures
%        |-- .CIupperbound - upper and lower (1 - alpha) confidence interval of 
%        |-- .CIlowerbound   |PDC|^2 or |DTF|^2 estimates
%        |-- .measure   - Type of Phi: 'pdc','gpdc','ipdc','dtf','dc' or 'idtf' 
%        |-- .SS        - Spectral density
%        |-- .Comment   - Comment
%      or
%        .{Phi,Pval,Threshold,CIupperbound,CIlowerbound,measure,SS}
%
%        d: structure variable with xplot.m (KS+LAB) compatible fields:
%        |-- .pdc2/dtf2 - |PDC|^2 / |DTF|^2  estimates
%        |-- .pvalues   - p-values associated to PDC2 / DTF2  estimates. 
%        |-- .th        - |PDC|^2 / |DTF|^2  threshold value with (1-avalue) 
%        |                significance level.
%        |-- .pdc2_th/dtf2_th - significant values above threshold level
%        |                otherwise NaN in freqs
%        |-- .{ci1,ci2} - upper and lower (1 - alpha) confidence interval of 
%        |                |PDC|^2 / |DTF|^2  estimates
%        |-- .metric    - metric used for |PDC|^2 or |DTF|^2  estimation:
%        |                'euc', 'diag' or 'info'
%        |-- .alpha     - significance level
%        |-- .p         - VAR model order
%        |-- .Comment   - comment
%        |-- .SS        - power spectra
%        +-- .coh2      - squared spectral coherence
%      or
%        .{{pdc2/dtf2},pvalues,th,{pdc2_th/dtf2_th},ci1,ci2,metric,alpha,p,
%   
%=========================================================================%
% The End
%=========================================================================%

%% Compatibility problem with Array Sizes on rdivide operation
% Reproduced from <https://www.mathworks.com/help/matlab/ref/rdivide.html>
%
%% Version History 
% Introduced before R2006a
% 

%% R2018b: sum function with 'all' argument
% "S = sum(A,'all') computes the sum of all elements of A. This syntax is valid
% for MATLAB® versions R2018b and later."
% (<https://www.mathworks.com/help/matlab/ref/sum.html#d124e1470580>

%% R2016b: Implicit expansion change affects arguments for operators
% Behavior changed in R2016b
% 
% "Starting in R2016b with the addition of implicit expansion, some combinations
% of arguments for basic operations that previously returned errors now produce
% results. For example, one previously could not add a row and a column vector,
% but those operands are now valid for addition. In other words, an expression
% like [1 2] + [1; 2] previously returned a size mismatch error, but now it
% executes.
% 
% If your code uses element-wise operators and relies on the errors that MATLAB
% previously returned for mismatched sizes, particularly within a try/catch
% block, then your code might no longer catch those errors.
% 
% For more information on the required input sizes for basic array operations,
% see Compatible Array Sizes for Basic Operations." 
% (<https://www.mathworks.com/help/matlab/ref/rdivide.html>)

%% R2020b: Implicit expansion change affects duration arrays
% Behavior changed in R2020b
% 
% "Starting in R2020b, rdivide supports implicit expansion when the arguments
% are duration arrays. Between R2020a and R2016b, implicit expansion was
% supported only for numeric data types."
% (<https://www.mathworks.com/help/matlab/ref/rdivide.html>)

%% if elsif else end template for rdivide operation. KS
%          if versionNumberMatlab < 9.2  % MATLAB before release 2017a
% 
%          elseif versionNumberMatlab < 9.9 % MATLAB before release 2020a
% 
%          else
%                         
%          end

%% help icdf
%  "icdf Inverse cumulative distribution function for a specified distribution.
%     . . .   
%     It is a generic function that accepts a distribution by name, and not
%     compatible with Octave v 6.4. It is faster to use a more specialized
%     function when possible, such as NORMINV for the normal distribution and,
%     in this function, CHI2INV is used for the Chisquare distribution.
%
%     Version History: Introduced before R2006a"
%     (https://www.mathworks.com/help/stats/prob.normaldistribution.icdf.html)

%% Logs:
%  2022.11.22 - Edting compatibility with MATLAB versions
