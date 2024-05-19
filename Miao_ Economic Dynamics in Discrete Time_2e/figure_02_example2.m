%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x'=f(x)=a*x*(1-x);
% transcritical bifurcation
% Example2 in chapter 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;

xmin=-1;

x=xmin:0.01:1;

a1=1;
y1=a1*x.*(1-x);
x1=0; %steady state
x21=(a1-1)/a1; %steady state

a2=2;
y2=a2*x.*(1-x);
x1=0;
x22=(a2-1)/a2;

a3=0.5;
y3=a3*x.*(1-x);
x1=0;
x23=(a3-1)/a3;



figure;
hold on;

z= plot(x,y1,x,y2,'--',x,y3,'-.','LineWidth',1.5)

hold on
plot(x,0*x);
plot(x,x);
xlabel('$x$','interpreter','latex','Fontsize',14);

ylabel('$f(x;a)$','interpreter','latex', 'Fontsize',14);


f =@ (x) a2*x.*(1-x);
cobweb(f,x22+0.2,4)
hold on
cobweb(f,x1+0.05,4)

f =@ (x) a3*x.*(1-x);
cobweb(f,x23+0.1,4)
cobweb(f,x1-0.2,4)

axis([xmin 1 xmin 1])
legend(z,'a=1','a=2','a=0.5','Location','SouthEast')
%print -depsc figure_02_example2.eps
