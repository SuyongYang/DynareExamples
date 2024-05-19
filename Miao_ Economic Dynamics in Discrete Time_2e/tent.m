%
% This program plots the tent map
%
close all;

x=0:0.0001:0.999;

f=@(x) 2*(0.5-abs(x-0.5));

figure
subplot(1,3,1);
plot(x,f(x),x,x,'--','LineWidth',1.5)
title('$f(x)$','interpreter','latex')
subplot(132)
plot(x,f(f(x)),x,x,'--','LineWidth',1.5)
title('$f^2(x)$','interpreter','latex')
subplot(133)
plot(x,f(f(f(x))),x,x,'--','LineWidth',1.5)
title({'$f^3(x)$'},'interpreter','latex')

print -depsc 9703_002_fig_009.eps