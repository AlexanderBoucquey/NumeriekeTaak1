function [lambda,v] = mmachten(A,v,maxit)

% function [lambda,v]= mmachten(A,v,maxit)
% 
% berekent met de methode van de machten de grootste eigenwaarde van de matrix A
%
% invoer
% A - matrix
% v - startvector
% maxit - maximum aantal iteraties
%
% uitvoer 
% lambda - de berekende eigenwaarde
% v - de berekende eigenvector

[n,m] = size(A);
if n~=m,
  disp('A is geen vierkante matrix')
  return
end
if n<2
  disp('A moet minstens dimensie 2 hebben')
  return
end

v = v/norm(v);

for it = 1:maxit
    w = A*v;
    v = w/norm(w);
    lambda = v'*A*v;
end
