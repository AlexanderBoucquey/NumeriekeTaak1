function [e,res]=qr_zonder(A)

% function [e]=qr_zonder(A)
%
% berekent met de QR-methode zonder shift een eigenwaarde van de matrix A
%
% invoer
% A - matrix
% 
% uitvoer
% e - de berekende eigenwaarde
% res - de normen van de residu's voor iedere iteratiestap

[n,m] = size(A);
if n~=m,
  disp('A is geen vierkante matrix')
  return
end
if n<2
  disp('A moet minstens dimensie 2 hebben')
  return
end

res = [];

while abs(A(n,n-1))>1.e-13
   res = [res abs(A(n,n-1))];
   [q,r]=qr(A);
   A = r*q;
end
res = [res abs(A(n,n-1))];
disp(sprintf('residu = %.1e', abs(A(n,n-1))))
e = A(n,n);
