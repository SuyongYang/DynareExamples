%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plots stable and unstable manifold and vector fields for the
% linear model in my second edition book
%              x(t+1)=-5/4*x(t)-3/4*y(t)
%              y(t+1)= -3/4*x(t)+5/4*y(t)
% Call the Matlab code: vectorfield.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;
clear all;

x=-1:0.1:1;

figure

z=plot(x,x/3,x,3*x,'--', x,x, '-kx', x, -x,'-mo','LineWidth',2.2)  


hold on


sys = @(t,x) [x(1)/4-3*x(2)/4; -3*x(1)/4+x(2)/4];
vectorfield(sys,-1:.2:1, -2:.2:2);
ylabel('y_t','FontSize',14);
xlabel('k_t','FontSize',15);
legend(z,'\Delta k_t=0','\Delta y_t=0','stable manifold','unstable manifold','Location','SouthEast')

print -depsc phase3.eps




