%% CMLSM
%        Multivariate autoregressive model least squares estimator.
%% Syntax
%       [npf,na,nef] = CMLSM(u,IP)
%
%% Input arguments
%         u   - vector of rows
%         IP  - order
%
%% Output arguments
%         npf - error covariance
%         na  - model
%         nef - residue

%% Code
function [npf,na,nef] = cmlsm(u,IP)

[m,n] = size(u);
[b,SU,nfe] = mlsmx(u,IP);
na = reshape(b,m,m,IP);
npf = SU*n; % see normalization

%==========================================================================
%
% 1998/01/30 (L.A.B.); 2000/04/11 (LAB reviewed)

function [b,SU,nfe] = mlsmx(Y,p)
[K,T] = size(Y);
Z = zmatrm(Y,p);
Gamma = Z*Z';
U1 = Gamma\Z; % Equivalent to inv(Gamma)*Z;

SU = (Y*Y'-Y*Z'*U1*Y');
SU = SU/(T-K*p-1);
b = kron(U1,eye(K))*reshape(Y,K*T,1);
nfe = reshape(reshape(Y,K*T,1)-kron(Z',eye(K))*b,K,T); 

%==========================================================================
% Computation of Z - data structure (no estimation of the mean)
%
% function Z = zmatrm(Y,p);
%
% input:  Y - data in row vectors 
%         p - model covariance order
%
% output: Z
%
% [1998.01.30]: L.A.B.
%
function Z = zmatrm(Y,p)
[K,T] = size(Y);
y1 = [zeros(K*p,1);reshape(flipud(Y),K*T,1)];
Z =  zeros(K*p,T);
for i = 0:T-1
   Z(:,i+1) = flipud(y1(1+K*i:K*i+K*p));
end
