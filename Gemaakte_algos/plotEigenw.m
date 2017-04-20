function plotEigenw(alpha)
m = 100;
A = sprand(m,m,0.5);
A = A + alpha*speye(m); A=A/norm(A,1);
b = rand(m,1);
[x, itx, res, res2] = GMRES(A,b);
plot(res);
set(gca, 'YScale', 'log');
pause;
plot(res2);
set(gca, 'YScale', 'log');
pause;
[B, C] = eigs(A);
figure;
hold on
r = 1;
th = 0:pi/50:2*pi;
xunit = r * cos(th);
yunit = r * sin(th);
h = plot(xunit, yunit);
plot(C, '*');
hold off




