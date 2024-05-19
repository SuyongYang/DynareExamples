% logbif.m - this MATLAB file simulates the 
% logistic difference equation
%  x(t+1)=a*x(t)*(1-x(t))
% and carries out a bifurcation analysis by varying a.
%
% 200 different values of a are used between the 
% ranges amin and amax set by the user. A bifurcation
% plot is drawn by showing the last 250 points of
% a sequence of n simulated points for each
% value of a in [amin, amax]. The initial conditions are 
% fixed at x0=0.3

close all
%clear all

amin=2.6;
amax=4;
x0=0.3;
n=900;

jmax=200;
t=zeros(jmax+1,1);
z=zeros(jmax+1,250);
del=(amax-amin)/jmax;
for j=1:jmax+1
x=zeros(n+1,1);
x(1)=x0;
t(j)=(j-1)*del+amin;
a=t(j);
 for i=1:n
  x(i+1)=a*x(i)*(1-x(i));
   if (i>750) 
   z(j,i-750)=x(i+1);
   end
 end
end
plot(t,z,'r.','MarkerSize',5)
xlabel('a','FontSize',12), 
ylabel('long-run values of  \itx_{t}','FontSize',12)
%title('Bifurcation diagram for the logistic difference equation')
print -depsc figure_01_bifur.eps