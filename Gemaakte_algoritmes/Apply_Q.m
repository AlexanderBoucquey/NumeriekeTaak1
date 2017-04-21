function y = Apply_Q(L,b)
sz = size(L);
m = sz(1);
n = sz(2);
for k= 1:n
    vk = L(k:n,k);
    b(k:m,1) = b(k:m,1) - 2* vk* (vk.'*b(k:m,1));
end
y=b;