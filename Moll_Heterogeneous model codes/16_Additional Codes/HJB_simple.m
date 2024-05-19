clear all; clc;

tic;

s = 2;
r = 0.045;
rho = 0.05;
w = .1;

I=500;
amin = -0.02;
amax = 1;
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);

maxit=10000;
crit = 10^(-6);

dVf = zeros(I,1);
dVb = zeros(I,1);
c = zeros(I,1);

%INITIAL GUESS
v0 = (w + r.*a).^(1-s)/(1-s)/rho;
v = v0;


for n=1:maxit
    V = v;
    % backward difference
    dVb(2:I) = (V(2:I)-V(1:I-1))/da;
    dVb(1) = (w + r.*amin).^(-s); %boundary condition
    
    I_concave = dVb > dVf; %indicator whether value function is concave (problems arise if this is not the case)
    
    c = dVb.^(-1/s);
    Vchange = c.^(1-s)/(1-s) + dVb.*(w + r.*a - c) - rho.*V;
       
    %% This is the update
    % the following CFL condition seems to work well in practice
    Delta = .9*da/max(w + r.*a);
    v = v + Delta*Vchange;
    
    dist(n) = max(abs(Vchange));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

% Graphs
set(gca,'FontSize',14)
plot(dist,'LineWidth',2)
grid
xlabel('Iteration')
ylabel('||V^{n+1} - V^n||')

Verr = c.^(1-s)/(1-s) + dVb.*(w + r.*a - c) - rho.*V;

set(gca,'FontSize',14)
plot(a,Verr,'LineWidth',2)
grid
xlabel('k')
ylabel('Error in HJB Equation')
xlim([amin amax])

adot = w + r.*a - c;

set(gca,'FontSize',12)
plot(a,V,'LineWidth',2)
grid
xlabel('a')
ylabel('V(a)')
xlim([amin amax])

set(gca,'FontSize',14)
plot(a,c,'LineWidth',2)
grid
xlabel('a')
ylabel('c(a)')
xlim([amin amax])

set(gca,'FontSize',14)
plot(a,adot,a,zeros(1,I),'--','LineWidth',2)
grid
xlabel('a')
ylabel('s(a)')
xlim([amin amax])
print -depsc HJB_stateconstraint_simple.eps

%Approximation at borrowing constraint

u1 = (w+r*amin)^(-s); u2 = -s*(w+r*amin)^(-s-1);
nu = sqrt(-2*(rho-r)*u1/u2);
s_approx = -nu*(a-amin).^(1/2);

set(gca,'FontSize',14)
h1 = plot(a,adot,a,s_approx,'-.',a,zeros(1,I),'--','LineWidth',2)
legend(h1,'s_1(a)','Approximate s_1(a)','Location','SouthWest')
grid
xlabel('a')
ylabel('s(a)')
xlim([amin amax])
%print -depsc HJB_stateconstraint_simple.eps