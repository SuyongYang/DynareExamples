%Optimized for speed by SeHyoun Ahn

clear all; clc; close all;

tic;

s = 2;
rho = 0.05;
z1 = .2;
z2 = 2*z1;
z = [z1,z2];
la1 = 1;
la2 = 1;
la = [la1,la2];
z_ave = (z1*la2+z2*la1)/(la1+la2);

Aprod = 0.3;
al = 1/3;
d = 0.05;


I=500;
amin = 0;
amax = 20;
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);

aa = [a,a];
zz = ones(I,1)*z;


maxit= 100;
crit = 10^(-6);
Delta = 1000;

dVf = zeros(I,2);
dVb = zeros(I,2);
c = zeros(I,2);

Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

Ir = 100;
rmin = -0.0499;
rmax = 0.049;
r_grid = linspace(rmin,rmax,Ir);

%INITIAL GUESS
r = r_grid(1);
KD = (al*Aprod/(r+d))^(1/(1-al))*z_ave;
w = (1-al)*Aprod*(KD/z_ave)^al;
% v0(:,1) = (z(1) + r.*a).^(1-s)/(1-s)/rho;
% v0(:,2) = (z(2) + r.*a).^(1-s)/(1-s)/rho;
v0(:,1) = (w*z(1) + max(r,0.01).*a).^(1-s)/(1-s)/rho;
v0(:,2) = (w*z(2) + max(r,0.01).*a).^(1-s)/(1-s)/rho;

for ir=1:Ir;

r = r_grid(ir);
KD(ir) = (al*Aprod/(r+d))^(1/(1-al))*z_ave;
w = (1-al)*Aprod*(KD(ir)/z_ave)^al;
w_r(ir)=w;


if ir>1
v0 = V_r(:,:,ir-1);
end

v = v0;

for n=1:maxit
    V = v;
    V_n(:,:,n)=V;
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (w*z + r.*amax).^(-s); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (w*z + r.*amin).^(-s); %state constraint boundary condition
    
    I_concave = dVb > dVf; %indicator whether value function is concave (problems arise if this is not the case)
    
    %consumption and savings with forward difference
    cf = dVf.^(-1/s);
    ssf = w*zz + r.*aa - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/s);
    ssb = w*zz + r.*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = w*zz + r.*aa;
    dV0 = c0.^(-s);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = ssf > 0; %positive drift --> forward difference
    Ib = ssb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    %make sure backward difference is used at amax
    %Ib(I,:) = 1; If(I,:) = 0;
    %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
    %already taken care of automatically
    
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = dV_Upwind.^(-1/s);
    u = c.^(1-s)/(1-s);
    
    %CONSTRUCT MATRIX
    X = - min(ssb,0)/da;
    Y = - max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;
    
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
b = zeros(2*I,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg = AT\b;
g_sum = gg'*ones(2*I,1)*da;
gg = gg./g_sum;

g = [gg(1:I),gg(I+1:2*I)];

check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;

g_r(:,:,ir) = g;
adot(:,:,ir) = w*zz + r.*aa - c;
V_r(:,:,ir) = V;
dV_r(:,:,ir) = dV_Upwind;
c_r(:,:,ir) = c;

S(ir) = g(:,1)'*a*da + g(:,2)'*a*da;
end

Smax = max(S);
amin1 = amin-0.02;
aaa = linspace(amin1,Smax,Ir);
rrr = linspace(-0.06,0.06,Ir);
KD = (al*Aprod./(max(rrr+d,0))).^(1/(1-al))*z_ave;

set(gca,'FontSize',14)
plot(S,r_grid,KD,rrr,zeros(1,Ir)+amin,rrr,'--',aaa,ones(1,Ir)*rho,'--',aaa,ones(1,Ir)*(-d),'--','LineWidth',2)
text(0.05,0.052,'$r = \rho$','FontSize',16,'interpreter','latex')
text(0.05,-0.054,'$r = -\delta$','FontSize',16,'interpreter','latex')
text(0.1,0.02,'$S(r)$','FontSize',16,'interpreter','latex')
text(0.29,0.02,'$r=F_K(K,1)-\delta$','FontSize',16,'interpreter','latex')
text(0.01,0.035,'$a=\underline{a}$','FontSize',16,'interpreter','latex')
ylabel('$r$','FontSize',16,'interpreter','latex')
xlabel('$K$','FontSize',16,'interpreter','latex')
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
ylim([-0.06 0.06])
xlim([amin1 0.6])
print -depsc aiyagari_asset_supply.eps