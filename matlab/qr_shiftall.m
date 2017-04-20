function [E,Res] = qr_shiftall(A)

% function [E,Res] = qr_shiftall(A)
%
% berekent met de QR-methode al de eigenwaarden van de matrix A
%
% invoer
% A - matrix
% 
% uitvoer
% E - de eigenwaarden van A
% Res - matrix van residu's (de residu's voor element A(i,i) staan op R(i,:))
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
Res = [];
E = [];
it = 1;
for i = n:-1:2
    while abs(A(i,i-1))>1.e-13
        mu = A(i,i);
        [q,r]=qr(A(1:i,1:i)-mu*eye(i));
        A(1:i,1:i) = r*q + mu*eye(i);
        Res(:,it)=abs(diag(A,-1));
        it = it +1;
    end
    E = [E A(i,i)];
end
E = [E A(1,1)];
