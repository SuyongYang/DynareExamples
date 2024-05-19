% Plot the logistic function
% y=f(x)=a*x*(1-x);
% and its third iterate.
%

close all;
clear all;

a=4;

x=0:0.01:1;

f=@(x) a*x.*(1-x);
figure
plot(x,f(x),'--',x,x,'-.',x,f(f(f(x))),'LineWidth',1.5)

xlabel('$x$','interpreter','latex','Fontsize',12);
legend({'$f(x)$','$x$','$f^3(x)$'},'interpreter','latex','Fontsize',12,'Location','SouthEast')

print -depsc 9703_002_fig_006.eps