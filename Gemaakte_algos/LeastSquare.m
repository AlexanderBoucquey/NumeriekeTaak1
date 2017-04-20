function [x, r] = LeastSquare(y, u, n)
sz = size(y);
m = sz(1);
A = zeros(m-n,n);
for k = 1:m-n
    A(k,:)= flipud(u(k:k+n-1));
end
b = y(n+1:m);
[Q, R] = Householder_explicit(A);
x = linsolve(R, Q.'*b);
r = norm(b-A*x);
% data = [x.' r];
% fileID = fopen('Output.txt','a');
% fprintf(fileID,'%6s','x');
% fprintf(fileID,'%8.4f',data.');
% fclose(fileID);

