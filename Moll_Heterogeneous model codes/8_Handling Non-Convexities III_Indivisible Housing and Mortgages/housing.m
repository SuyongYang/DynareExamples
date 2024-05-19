clear all; clc; close all;

tic;

s = 2;
r = 0.035;
al = 1/3;
eta = 0.2;
p = 1;
phi = 2;
hmin = 2.3;
p*hmin/phi
rho = 0.05;
z1 = .1;
z2 = .2;
z2 = .135;
z = [z1,z2];
la1 = 0.5;
la2 = 0.5;
la = [la1,la2];


I=500;
amin = 0;
amax = 3;
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);

aa = [a,a];
zz = ones(I,1)*z;


maxit= 120;
crit = 10^(-6);
Delta = 1000;

dVf = zeros(I,2);
dVb = zeros(I,2);
c = zeros(I,2);

Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

h = min((al*eta/(r*p))^(1/(1-al)) + hmin,phi*aa/p);
h = h.*(h>=hmin);
f = eta*(max(h-hmin,0)).^al - r*p*h;

% h = min((al*eta/(r*q))^(1/(1-al)),phi*aa/q);
% h = h.*(h>=hmin);
% f = eta*h.^al - r*q*h;

%INITIAL GUESS
%v0 = (zz + f + r.*aa).^(1-s)/(1-s)/rho;
v0 = (zz + r*aa).^(1-s)/(1-s)/rho;

v = v0;


for n=1:maxit
    V = v;
    V_n(:,:,n)=V;
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (z + f(I,:) + r.*amax).^(-s); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (z + f(1,:) + r.*amin).^(-s); %state constraint boundary condition
       
    %consumption and savings with forward difference
    cf = (max(dVf,10^(-10))).^(-1/s);
    ssf = zz + f + r.*aa - cf;
    Hf = cf.^(1-s)/(1-s) + dVf.*ssf;
    %consumption and savings with backward difference
    cb = (max(dVb,10^(-10))).^(-1/s);
    ssb = zz + f + r.*aa - cb;
    Hb = cb.^(1-s)/(1-s) + dVb.*ssb;
    %consumption and derivative of value function at steady state
    c0 = zz + f + r.*aa;
    
%     % dV_upwind makes a choice of forward or backward differences based on    
    Ineither = (1-(ssf>0)) .* (1-(ssb<0));
    Iunique = (ssb<0).*(1-(ssf>0)) + (1-(ssb<0)).*(ssf>0);
    Iboth = (ssb<0).*(ssf>0);
    Ib = Iunique.*(ssb<0) + Iboth.*(Hb>=Hf);
    If = Iunique.*(ssf>0) + Iboth.*(Hf>=Hb);
    I0 = Ineither;
    
    c = cf.*If + cb.*Ib + c0.*I0;
    u = c.^(1-s)/(1-s);
    
    %CONSTRUCT MATRIX
    X = - Ib.*ssb/da;
    Y = - If.*ssf/da + Ib.*ssb/da;
    Z = If.*ssf/da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);

    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    B = (1/Delta + rho)*speye(2*I) - A;
    
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%
AT = A';

%INITIAL DISTRIBUTION 1:
gg0 = ones(2*I,1);

%INITIAL DISTRIBUTION 2
%j = I/2; gg0 = [zeros(j,1);ones(I-j,1);zeros(j,1);ones(I-j,1)];

%INITIAL DISTRIBUTION 3
%j = 10; gg0 = [ones(j,1);zeros(I-j,1);ones(j,1);zeros(I-j,1)];

g_sum = gg0'*ones(2*I,1)*da;
gg0 = gg0./g_sum;

gg{1}=gg0;
N=1000; dt=10;
for n=1:N
    %Implicit method in Updating Distribution.
    gg{n+1}= (speye(2*I) - AT*dt)\gg{n};
    g_dist(n)=max(abs(gg{n+1}-gg{n}));
end

g = [gg{N}(1:I),gg{N}(I+1:2*I)];


adot = zz + f + r.*aa - c;

astar = p*hmin/phi;
[obj, index] = min(abs(astar-a));


amax1 = amax;
amin1 = -0.1;

figure(1)
h1 = plot(a,adot(:,1),'b',a,adot(:,2),'r-.','LineWidth',2)
legend(h1,{'$s_1(a)$','$s_2(a)$'},'FontSize',20,'interpreter','latex')
hold on
plot(linspace(amin1,amax,I),zeros(I,1),'k--')
set(gca,'FontSize',16)
set(gca, 'XTick', []);
set(gca, 'YTick', []);
line([astar astar], [-0.08 0.1],'Color','Black','LineStyle','--')
line([amin amin], [-0.08 0.1],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex','FontSize',20)
ylabel('Saving, $s_j(a)$','interpreter','latex','FontSize',20)
xlim([amin1 amax1])
print -depsc saving_housing.eps
hold off

c = c - eta*(max(h-hmin,0)).^al;

figure(2)
h1 = plot(a,c(:,1),'b',a,c(:,2),'r-.','LineWidth',2)
set(gca,'FontSize',16)
set(gca, 'XTick', []);
set(gca, 'YTick', []);
legend(h1,{'$c_1(a)$','$c_2(a)$'},'FontSize',20,'interpreter','latex')
line([astar astar], [0 0.17],'Color','Black','LineStyle','--')
line([amin amin], [0 0.17],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex','FontSize',20)
ylabel('Consumption, $c_j(a)$','interpreter','latex','FontSize',20)
xlim([amin1 amax1])
ylim([0 0.17])
print -depsc consumption_housing.eps

figure(3)
h1 = plot(a,h(:,1),'b',a,h(:,2),'r--','LineWidth',2)
set(gca,'FontSize',16)
set(gca, 'XTick', []);
set(gca, 'YTick', []);
line([astar astar], [0 5],'Color','Black','LineStyle','--')
line([amin amin], [0 5],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex','FontSize',20)
ylabel('Housing, $h_j(a)$','interpreter','latex','FontSize',20)
xlim([amin1 amax1])
print -depsc housing_housing.eps

figure(4)
h1 = plot(a,f(:,1),'b',a,f(:,2),'r--','LineWidth',2)
hold on
plot(linspace(amin1,amax,I),zeros(I,1),'k--')
set(gca,'FontSize',16)
set(gca, 'XTick', []);
set(gca, 'YTick', []);
line([amin amin], [-0.04 0.12],'Color','Black','LineStyle','--')
line([astar astar], [-0.04 0.12],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex','FontSize',20)
ylabel('Pecuniary Benefit from Housing, $f(a)$','interpreter','latex','FontSize',20)
xlim([amin1 amax1])
print -depsc payoff_housing.eps
hold off

figure(5)
h1 = plot(a,v(:,1),'b',a,v(:,2),'r-.','LineWidth',2)
set(gca,'FontSize',16)
set(gca, 'XTick', []);
set(gca, 'YTick', []);
legend(h1,{'$v_1(a)$','$v_2(a)$'},'FontSize',20,'interpreter','latex','Location','SouthEast')
line([astar astar], [-180 -60],'Color','Black','LineStyle','--')
line([amin amin], [-180 -60],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex','FontSize',20)
ylabel('Value Function, $v_j(a)$','interpreter','latex','FontSize',20)
xlim([amin1 amax1])
ylim([-180 -60])
print -depsc value_housing.eps

% plot(a,v(:,1),'b',a,v(:,2),'r','LineWidth',2)
% xlim([0 1])

figure(6)
h1 = plot(a,g(:,1),'b',a,g(:,2),'r-.','LineWidth',2)
set(gca,'FontSize',16)
set(gca, 'XTick', []);
set(gca, 'YTick', []);
legend(h1,{'$g_1(a)$','$g_2(a)$'},'FontSize',20,'interpreter','latex')
line([astar astar], [0 3.5],'Color','Black','LineStyle','--')
line([amin amin], [0 3.5],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex','FontSize',20)
ylabel('Densities, $g_j(a)$','interpreter','latex','FontSize',20)
xlim([amin1 amax1])
ylim([0 3.5])
print -depsc densities_housing.eps



% G(:,1)=cumsum(g(:,1)*da);G(:,2)=cumsum(g(:,2)*da);
% 
% figure(6)
% h1 = plot(a,G(:,1),'b',a,G(:,2),'r','LineWidth',2)
% set(gca,'FontSize',16)
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% legend(h1,'G_1(a)','G_2(a)','Location','SouthEast')
% line([astar astar], [0 0.5],'Color','Black','LineStyle','--')
% line([amin amin], [0 0.5],'Color','Black','LineStyle','--')
% xlabel('Wealth, $a$','interpreter','latex')
% ylabel('Cumulative Distribution Functions, $G_j(a)$','interpreter','latex')
% xlim([amin1 amax1])
% %print -depsc CDFs_housing.eps


%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%
AT = A';

%INITIAL DISTRIBUTION 1:
%gg0 = ones(2*I,1);

%INITIAL DISTRIBUTION 2
j = I/2; gg0 = [zeros(j,1);ones(I-j,1);zeros(j,1);ones(I-j,1)];

%INITIAL DISTRIBUTION 3
%j = 10; gg0 = [ones(j,1);zeros(I-j,1);ones(j,1);zeros(I-j,1)];

g_sum = gg0'*ones(2*I,1)*da;
gg0 = gg0./g_sum;

gg{1}=gg0;
N=1000; dt=10;
for n=1:N
    %Implicit method in Updating Distribution.
    gg{n+1}= (speye(2*I) - AT*dt)\gg{n};
    g_dist(n)=max(abs(gg{n+1}-gg{n}));
end

g = [gg{N}(1:I),gg{N}(I+1:2*I)];


figure(7)
h1 = plot(a,g(:,1),'b',a,g(:,2),'r-.','LineWidth',2)
set(gca,'FontSize',16)
set(gca, 'XTick', []);
set(gca, 'YTick', []);
legend(h1,{'$g_1(a)$','$g_2(a)$'},'FontSize',20,'interpreter','latex')
line([astar astar], [0 3.5],'Color','Black','LineStyle','--')
line([amin amin], [0 3.5],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex','FontSize',20)
ylabel('Densities, $g_j(a)$','interpreter','latex','FontSize',20)
xlim([amin1 amax1])
ylim([0 3.5])
print -depsc densities_housing2.eps
