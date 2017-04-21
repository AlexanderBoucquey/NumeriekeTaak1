function [E] = bisection_eigenvalue(A, a, b, tol)
[m, n] = size(A);
E_old =[];
for i = 1:m
    p = poly(A(1:i, 1:i));
    [k,l] = size(E_old);
    E = [];
    a_new = a;
    for j = 1:l+1
        if (l == 0)
            x = bisection(p,a,b,tol);
            E = [E x];
        elseif (j == l+1)               
            x = bisection(p,a_new,b,tol);  
%             disp(x);
            E = [E x]; 
        else
            x = bisection(p,a_new,E_old(1,j),tol);
            a_new = E_old(1,j);
            E = [E x];
        end    
    end    
    E_old = E;
    disp(E_old);
end
