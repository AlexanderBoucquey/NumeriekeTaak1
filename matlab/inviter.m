function [lambda, v] = inviter(A, mu, v, maxit);

% function [lambda, v] = inviter(A, mu, v, maxit);
% 
% inverse iteratie
%
% invoer
% A   - symmetrische matrix
% mu  - benadering voor gezochte eigenwaarde
% v   - startvector
% maxit - maximum aantal iteraties
%
% uitvoer
% lambda - eigenwaarde
% v      - eigenvector


[n,m] = size(A);
if n~=m,
  disp('A is geen vierkante matrix')
  return
end
if n<2
  disp('A moet minstens dimensie 2 hebben')
  return
end

if length(v)~=n
    disp('v moet evenveel rijen als A hebben')
    return
end

v = v / norm(v);
lambda = mu;

for it = 1:maxit
  v = (A - mu * speye(n)) \ v;
  v = v / norm(v);
  lambda = v'*A*v ;
end

