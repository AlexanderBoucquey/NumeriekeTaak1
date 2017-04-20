function [x, itx, res, res2] = GMRES(A, b)
Qn(:,1) = b/norm(b);
r = 1;
itx = [];
n = 1;
[i, k] = size(A);
x = zeros(k, 1);
res = [];
res2 = [];
z = A\b;
while r > 10^-12
  e1 = zeros(n+1,1);
  e1(1,1) = 1;
  v = A*Qn(:,n);
  for j = 1:n
    H(j,n) = Qn(:,j)'*v;
    v = v - H(j,n)*Qn(:,j);
  end
  H(n+1,n) = norm(v);
  if H(n+1,n) <= 0, 
      break;
  end
  Qn(:,n+1) = v/H(n+1,n);
  [Q, R] = qr(H,0);
  y = linsolve(R, Q.'*(norm(b)*e1));
  x(:,1) = Qn(:,1:n)*y;
  itx = [itx x];
  r = norm(b-A*x); 
  res = [res r];
  res2 = [res2 norm(x-z)];
  n = n +1;
end;

