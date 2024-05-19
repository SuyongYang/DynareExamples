clear all; clc; close all;

rho = 0.05;

I= 1000;
xmin = 0.1; %Can't be zero, have mu_bar negative, and have sigma(xmin) = 0.  Otherwise absorbing state
xmax=1;
x = linspace(xmin,xmax,I)';
dx = x(2)-x(1);
dx2 = dx^2;

u = x.^(0.5);
mu_bar = -0.01;
mu = ones(I,1).*mu_bar;

sig_bar = 0.01;
sig = x.*sig_bar;
sig2 = sig.^2;

S = 10*ones(I,1);
%S = 10*exp(0.05*x);


%CONSTRUCT MATRIX
X = - min(mu,0)/dx + sig2/(2*dx2);
Y = - max(mu,0)/dx + min(mu,0)/dx - sig2/dx2;
Z =  max(mu,0)/dx + sig2/(2*dx2);

A =spdiags(Y,0,I,I)+spdiags(X(2:I),-1,I,I)+spdiags([0;Z(1:I-1)],1,I,I);

%Manually adjust the boundary values at the corners.
A(1,1) = Y(1) + X(1); %Reflecting barrier
A(I,I) = Y(I) + Z(I); %Reflecting barrier


%Solve option problem as LCP
B = rho*speye(I) - A;
q = -u + B*S; 

%TRY DIFFERENT ALGORITHMS FOR SOLVING LCP
tic
disp('Solving LCP')

%Test 1: Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
z0 = zeros(I,1); l = zeros(I,1); u = Inf*ones(I,1);
z = LCP(B,q,l,u,z0,1);

%Test 2: Andreas Almqvist pivoting (Lemke) LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/41485
%[w,z,retcode] = LCPSolve(B,q);

LCP_error = max(abs(z.*(B*z + q)));
if LCP_error > 10^(-6)
    disp('LCP not solved')
    break
end
    
V_LCP = z+S; %calculate value function
toc

error = z.*(B*z + q);
plot(x,error)

plot(x,V_LCP,x,S,'--','LineWidth',2)
set(gca,'FontSize',16)
legend('v(x)','S(x)','Location','NorthWest')
xlabel('x')
%print -depsc option_simple.eps
