function [e,res]=qr_shiftrayleigh(A)

% function [e,res]=qr_shiftrayleigh(A)
%
% berekent met de QR-methode met Rayleigh quotient shift een eigenwaarde van de matrix A
%
% invoer
% A - matrix
% 
% uitvoer
% e - de berekende eigenwaarde
% res - de normen van de residu's voor iedere iteratiestap
%
% De gebruikte methode is de QR-methode met shift. Als shift wordt het
% (n,n)-element van A gebruikt.

[n,m] = size(A);
if n~=m,
  disp('A is geen vierkante matrix')
  return
end
if n<2
  disp('A moet minstens dimensie 2 hebben')
  return
end

A = hess(A);
res = [];

while abs(A(n,n-1))>1.e-13
   res = [res abs(A(n,n-1))];
   [q,r]=qr(A-A(n,n)*eye(n));
   A = r*q + A(n,n)*eye(n);
end
res = [res abs(A(n,n-1))];
disp(sprintf('residu = %.1e', abs(A(n,n-1))));
e = A(n,n);
