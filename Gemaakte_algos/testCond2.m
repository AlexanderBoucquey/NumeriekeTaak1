function [nx, nr] = testCond2(m, k)
A = rand(m);
b = rand(m,1);
[U, S, V] = svd(A);
S(1,1) = S(m, m)*k;
if k == 1
    S = eye(m);
end
A = U*S*V.';
disp(cond(A));
z= linsolve(A,b);
[L,R]= Householder_implicit(A);y = Apply_Q(L,b); x = linsolve(R,y);
dx = x-z;
r = A*x - b;
nx = norm(dx)/norm(x);
nr = norm(r)/norm(b);