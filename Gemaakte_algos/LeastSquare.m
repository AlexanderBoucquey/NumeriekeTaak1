function [x, r] = LeastSquare(y, u, u_validate, y_validate)
y_fault = [];
for n = 1:30
sz = size(y);
m = sz(1);
A = zeros(m-n,n);
for k = 1:m-n
    A(k,:)= flipud(u(k:k+n-1));
end
b = y(n+1:m);
[Q, R] = Householder_explicit(A);
x = linsolve(R, Q.'*b);

sz = size(y_validate);
m = sz(1);
disp(m);
A = zeros(m-n,n);
for k = 1:m-n
    A(k,:)= flipud(u_validate(k:k+n-1));
end
disp(size(A*x));
disp(size(x));
disp(size(y_validate(n+1:200,:)));
disp(norm(y_validate(n+1:200,:) - (A*x)));
y_fault = [y_fault norm(y_validate(n+1:200,:) - A*x)];

end

semilogy(y_fault);
% data = [x.' r];
% fileID = fopen('Output.txt','a');
% fprintf(fileID,'%6s','x');
% fprintf(fileID,'%8.4f',data.');
% fclose(fileID);

