function [H, Q, Z] = arnoldi(A, b, maxit)

% function [H, Q] = arnoldi(A, b, maxit)
%
% Arnoldi iteratie
% 
% invoer:
% A     - ijle matrix
% b     - startvector
% maxit - aantal iteraties
%
% uitvoer:
% H     - Hessenberg matrix
% Q     - orthogonale matrix


Q(:,1) = b/norm(b);
Z = zeros(4,maxit);
for n=1:maxit
  v = A*Q(:,n);
  for j = 1:n
    H(j,n) = Q(:,j)'*v;
    v = v - H(j,n)*Q(:,j);
  end
  d = eigs(H,min(n,4));
  H(n+1,n) = norm(v);
  if H(n+1,n) <= 0, 
      break;
  end
  Q(:,n+1) = v/H(n+1,n);
  Z(1: min(n,4),n) = d; 
end;

