function [Q, R] = Householder_explicit(A)
sz = size(A);
m = sz(1);
n = sz(2);
Q = eye(m);
for k = 1:n
    e1 = zeros(m-k+1,1);
    e1(1,1) = 1;
    x = A(k:m,k);
    vk = sign(x(1))*norm(x,2)*e1 + x;  
    F = eye(m-k+1) - 2*(vk*vk.')/(vk.'*vk);
    P = zeros(m);
    P(1:k-1,1:k-1) = eye(k-1);
    P(k:m,k:m) = F;
    A = P*A;
    Q = Q*P;    
end
R = A;