%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x'=f(x)=x^2+a;
% period doubling bifurcation
% Example3 in chapter 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
xmin=-1.2;
xmax=0.4;
x=xmin:0.01:xmax;
%steady state


a1=-3/4;
y1=x.^2+a1;
x21= 0.5+sqrt(1-4*a1)/2;
x11=0.5-sqrt(1-4*a1)/2;


a2=-3/4+0.1;
y2=x.^2+a2;
x22= 0.5+sqrt(1-4*a2)/2;
x12=0.5-sqrt(1-4*a2)/2;

a3=-3/4-0.1;
y3=x.^2+a3;
x23= 0.5+sqrt(1-4*a3)/2;
x13=0.5-sqrt(1-4*a3)/2;

a4=-1;
y4=x.^2+a4;
x24= 0.5+sqrt(1-4*a4)/2;
x14=0.5-sqrt(1-4*a4)/2;


figure;

z= plot(x,y2,x,y1,'--',x,y3,'-.',x,y4,':','LineWidth',1.5)

hold on
plot(x,0*x);
plot(x,x)
xlabel('$x$','interpreter','latex','Fontsize',14);

ylabel('$f(x;a)$','interpreter','latex', 'Fontsize',14);

n=7;

f =@ (x) x.^2+a2;
cobweb(f,x12+0.05,n)
hold on

f =@ (x) x.^2+a3;
cobweb(f,x13-0.05,n)

f =@ (x) x.^2+a4;
cobweb(f,x14-0.3,10)

cycle1=(-1-sqrt(-3-4*a4))/2;
cycle2=(-1+sqrt(-3-4*a4))/2;

axis([xmin xmax xmin xmax])
legend(z, 'a=-0.65','a=-0.75','a=-0.85','a=-1','Location','SouthEast')
%print -depsc figure_02_example3.eps
