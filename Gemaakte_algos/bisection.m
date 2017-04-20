function x = bisection(p, a, b, tol)
x = [];
while (abs(a-b) >= tol)
        fa = polyval(p, a);
        fb = polyval(p, b);    
        if (fa*fb > 0)
            disp('There is no root in the interval'); 
            x = [];
            break
        else 
            x = (a+b)/2;
            fx = polyval(p, x);           
            if (fx == 0) || (abs(fx) < tol)
                a = x;
                b = x;
            elseif (fa * fx < 0)
                b = x;
            else 
                a = x;
            end
        end
end