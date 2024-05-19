%
% This program plots the cobweb diagram for the logisitic map 
% f(x) = a*x(1-x) for a = 1.
% Call cobweb.m
%

close all
x=-0.1:0.01:1.1;

a=2;
x0=0.25;
f=@(x) a*x.*(1-x);
y1=f(x);
y2=x;
y3=0*x;
axes('FontSize',14);
plot(x,y1,'LineWidth',1.5)

hold on 
plot(x,y2,x,y3);

ylabel({'$x_{t+1}$'},'interpreter','latex','Fontsize',12)
xlabel({'$x_{t}$'},'interpreter','latex','Fontsize',12)
axis([-0.1 1.1 y1(1) y2(end)]);

cobweb(f,x0,7)

print -depsc 9703_001_fig_005.eps
