function plotArnoldi(A,Z)
[m,n] = size(Z);
[B, C] = eigs(A,4);
C = real(diag(C));
x = [1:150];
figure
plot(x, real(C(4,1)-Z(4,:)), '*')
set(gca,'YScale','log')
for i = 1:n
%     plot(C,'*')
%     hold on
%     plot(real(Z(:,i)),'*')
%     axis equal
%     hold off
%     pause;
%     plot(C(1,1)-Z(1,n))
%     
%     hold on
%     hold off
%     pause;

end

