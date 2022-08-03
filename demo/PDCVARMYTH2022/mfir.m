function  s = mfir(vv,b)
[m,m,q] = size(b);
[m,n,r] = size(vv);
s = reshape(vv,m,n*r);
for i = 1:q
   y = inzero(vv);
   s = s+b(:,:,i)*y;
   vv = reshape(y,m,n,r);
end
s = reshape(s,m,n,r);
