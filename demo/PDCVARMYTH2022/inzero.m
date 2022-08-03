function y = inzero(A)
[m,n,r] = size(A);
B = zeros(size(A));
B(:,2:n,:) = A(:,1:n-1,:);
y = reshape(B,m,n*r);