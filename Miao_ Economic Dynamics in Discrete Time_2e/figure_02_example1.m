%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x'=f(x)=x^2+a;
% saddle-node bifurcation
% Example1 in chapter 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;

x=0:0.01:1;

a1=0.5;
y1=x.^2+a1;

a2=1/4;
y2=x.^2+a2;

a3=1/7;

y3=x.^2+a3;

%figure;

z= plot(x,y1, '--',x,y2, x,y3, '-.','LineWidth',1.5)

hold on

plot(x,x);

xlabel('$x$','interpreter','latex','Fontsize',14);

ylabel('$f(x;a)$','interpreter','latex', 'Fontsize',14);
hold on


f =@ (x) x.^2+a3;
cobweb(f,0.05,4)
hold on


cobweb(f,0.3,4)

cobweb(f,0.8,4)

cobweb(f,0.85,4)

axis([0 1 0 1])
legend(z, 'a=1/2','a=1/4','a=1/7', 'Location','SouthEast')

%print -depsc C:\Users\4321`\documents\figure_02_example1.eps
