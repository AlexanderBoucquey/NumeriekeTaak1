function [L, R] = Householder_implicit(A)
sz = size(A);
m = sz(1);
n = sz(2);
L = zeros(m,n);
for k = 1:n
    e1 = zeros(m-k+1,1);
    e1(1,1) = 1;
    x = A(k:m,k);
    vk = sign(x(1))*norm(x,2)*e1 + x;    
    vk = vk / (norm(vk,2));
    A(k:m,k:n) = A(k:m,k:n) - 2* vk* (vk.'*A(k:m,k:n));
    L(k:m,k) = vk;
end
 R = A;