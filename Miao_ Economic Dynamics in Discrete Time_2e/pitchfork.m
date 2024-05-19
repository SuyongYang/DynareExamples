%Pitchfork bifurcation
% Plot map f(x)=ax-x^3
close all

x=-1:0.01:1;

a=1;

y1=a*x-x.^3;
y2=0.5*x-x.^3;
y3=1.5*x-x.^3;

figure
plot(x,y1,x,y2,'--',x,y3,'-.',x,x,':','LineWidth',1.5)
legend({'$a=1$','$a=0.5$','$a=1.5$','$y=x$'},'interpreter','latex','FontSize',12,'Location','SouthEast')

print -depsc 9703_002_fig_008.eps