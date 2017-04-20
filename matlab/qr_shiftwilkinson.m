function [e,res]=qr_shiftwilkinson(A)

% function [e,res]=qr_shiftwilkinson(A)
%
% berekent met de QR-methode een eigenwaarde van de matrix A
%
% invoer
% A - matrix
% 
% uitvoer
% e - de berekende eigenwaarde
% res - de normen van de residu's voor iedere iteratiestap
%
% De gebruikte methode is de QR-methode met shift. Als shift wordt het
% Wilkinson shift gebruikt.


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
   delta = (A(n-1,n-1)-A(n,n))/2;
   if(delta==0)
       signdelta = 1;
   else
       signdelta = sign(delta);
   end
   mu = A(n,n)-signdelta*A(n,n-1)^2/(abs(delta)+sqrt(delta^2+A(n,n-1)^2));
   [q,r]=qr(A-mu*eye(n));
   A = r*q + mu*eye(n);
end
res = [res abs(A(n,n-1))];
disp(sprintf('residu = %.1e', abs(A(n,n-1))));
e = A(n,n);
