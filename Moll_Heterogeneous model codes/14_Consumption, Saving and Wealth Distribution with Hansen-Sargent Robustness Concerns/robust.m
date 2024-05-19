clear all; close all; clc;
tic;
ga = 2; %CRRA utility with parameter gamma
r = 0.02; %interest rate
rho = 0.05; %discount rate

%ORNSTEIN-UHLENBECK PROCESS dlog(z) = -the*log(z)dt + sig2*dW
%STATIONARY DISTRIBUTION IS log(z) ~ N(0,Var) WHERE Var = sig2/(2*the)
Var = 0.07;
zmean = exp(Var/2); %MEAN OF LOG-NORMAL DISTRIBUTION N(0,Var)
Corr = 0.9;
nu = -log(Corr);
sig2 = 2*nu*Var;

theta = 0.5; %parameter

I=100;
amin = -0.02; %borrowing constraint
amax = 4;
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);

J=40;
zmin = zmean*0.8;
zmax = zmean*1.2;
z = linspace(zmin,zmax,J);
dz = (zmax-zmin)/(J-1);
dz2 = dz^2;

%plot(z,lognpdf(z,0,Var)) %PLOT STATIONARY DISTRIBUTION AND CHECK THAT CHOICE OF GRID DOESN'T CUT OFF TOO MUCH OF TAILS

aa = a*ones(1,J);
zz = ones(I,1)*z;

mu = (-nu*log(z) + sig2/2).*z; %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*z.^2; %VARIANCE (FROM ITO'S LEMMA)

s2z = ones(I,1)*s2;
sigz = sqrt(s2z);
muz = ones(I,1)*mu;

maxit= 100;
crit = 10^(-6);
Delta = 1000;

Vaf = zeros(I,J);
Vab = zeros(I,J);
Vzf = zeros(I,J);
Vzb = zeros(I,J);
Vzz = zeros(I,J);
c = zeros(I,J);

%INITIAL GUESS
v0 = (zz + r.*aa).^(1-ga)/(1-ga)/rho;
v = v0;

for n=1:maxit
    V = v;
    % forward difference
    Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    Vaf(I,:) = (z + r.*amax).^(-ga); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    Vab(1,:) = (z + r.*amin).^(-ga); %state constraint boundary condition
    
    I_concave = Vab > Vaf; %indicator whether value function is concave (problems arise if this is not the case)   
        
    %consumption and savings with forward difference
    cf = max(Vaf,10^(-8)).^(-1/ga);
    sf = zz + r.*aa - cf;
    %consumption and savings with backward difference
    cb = max(Vab,10^(-8)).^(-1/ga);
    sb = zz + r.*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = zz + r.*aa;
    
    % Upwind method makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = sf > 0; %positive drift --> forward difference
    Ib = sb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    
    c =  cf.*If + cb.*Ib + c0.*I0;
    
    %h
    Vzf(:,1:J-1) = (V(:,2:J) - V(:,1:J-1))/dz;
    Vzf(:,J) = theta.*muz(:,J)./s2z(:,J);
    Vzb(:,2:J) = (V(:,2:J) - V(:,1:J-1))/dz;
    Vzb(:,1) = theta.*muz(:,1)./s2z(:,1);
    
    hf = - Vzf.*sigz/theta;
    hb = - Vzb.*sigz/theta;
    h0 = - muz./max(sigz,10^(-8)); %h such that mu(z) + h*sigma(z) = 0
    
    Izf = (muz + sigz.*hf > 0);
    Izb = (muz + sigz.*hb < 0);
    Iz0 = 1 - Izf - Izb;
    
    h = hf.*Izf + hb.*Izb + h0.*Iz0;
    
    %utility function
    u = c.^(1-ga)/(1-ga) + theta/2*h.^2;
    
    
    %CONSTRUCT MATRIX A
    X = - min(sb,0)/da;
    Y = - max(sf,0)/da + min(sb,0)/da;
    Z = max(sf,0)/da;
        
    updiag=0; %This is needed because of the peculiarity of spdiags.
    for j=1:J
        updiag=[updiag;Z(1:I-1,j);0];
    end
    
    centdiag=reshape(Y,I*J,1);
    
    lowdiag=X(2:I,1);
    for j=2:J
        lowdiag=[lowdiag;0;X(2:I,j)];
    end
    
    AA=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);
    
    %CONSTRUCT MATRIX Bswitch SUMMARIZING EVOLUTION OF z
    chi_z =  - min(muz + hb.*sigz,0)/dz + s2z/(2*dz2);
    yy_z =  min(muz + hb.*sigz,0)/dz - max(muz + hf.*sigz,0)/dz - s2z/dz2;
    zeta_z = max(muz + hf.*sigz,0)/dz + s2z/(2*dz2);
       
    %This will be the upperdiagonal of the B_switch
    updiag = reshape(zeta_z,I*J,1);
    updiag = [zeros(I,1); updiag];

    %This will be the center diagonal of the B_switch    
    centdiag                = reshape(yy_z,I*J,1);
    centdiag(1:I)           = chi_z(:,1) + yy_z(:,1);
    centdiag(I*(J-1)+1:I*J) = yy_z(:,J) + zeta_z(:,J);    
    
    %This will be the lower diagonal of the B_switch
    lowdiag = reshape(chi_z(:,2:J),I*(J-1),1);

    %Add up the upper, center, and lower diagonal into a sparse matrix
    Bswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);  

    %Overall transition matrix for (a,z)
    A = AA + Bswitch;
    B = (1/Delta + rho)*speye(I*J) - A;
    
    u_stacked = reshape(u,I*J,1);
    V_stacked = reshape(V,I*J,1);
    
    b = u_stacked + V_stacked/Delta;

    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = reshape(V_stacked,I,J);
    
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

c_robust = c;
ss_robust = zz + r.*aa - c;
h_robust = h;
drift_robust = muz+h_robust.*sigz;

clear c h

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%NEED TO RECONSTRUCT A-matrix to make sure physical rather than perceived
%dynamics

%CONSTRUCT MATRIX Bswitch SUMMARIZING EVOLUTION OF z
chi =  - min(mu,0)/dz + s2/(2*dz2);
yy =  min(mu,0)/dz - max(mu,0)/dz - s2/dz2;
zeta = max(mu,0)/dz + s2/(2*dz2);

%This will be the upperdiagonal of the B_switch
updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
for j=1:J
    updiag=[updiag;repmat(zeta(j),I,1)];
end

%This will be the center diagonal of the B_switch
centdiag=repmat(chi(1)+yy(1),I,1);
for j=2:J-1
    centdiag=[centdiag;repmat(yy(j),I,1)];
end
centdiag=[centdiag;repmat(yy(J)+zeta(J),I,1)];

%This will be the lower diagonal of the B_switch
lowdiag=repmat(chi(2),I,1);
for j=3:J
    lowdiag=[lowdiag;repmat(chi(j),I,1)];
end

%Add up the upper, center, and lower diagonal into a sparse matrix
Bswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);

%CONSTRUCT MATRIX A
X = - min(sb,0)/da;
Y = - max(sf,0)/da + min(sb,0)/da;
Z = max(sf,0)/da;

updiag=0; %This is needed because of the peculiarity of spdiags.
for j=1:J
    updiag=[updiag;Z(1:I-1,j);0];
end

centdiag=reshape(Y,I*J,1);

lowdiag=X(2:I,1);
for j=2:J
    lowdiag=[lowdiag;0;X(2:I,j)];
end

AA=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);

A = AA + Bswitch;

%Solve Fokker-Planck Equation
AT = A';
b = zeros(I*J,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,I*J-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg = AT\b;
g_sum = gg'*ones(I*J,1)*da*dz;
gg = gg./g_sum;

g_robust = reshape(gg,I,J);

g_a_robust = sum(g_robust*dz,2);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOR COMPARISON: PROBLEM WITHOUT ROBUSTNESS CONCERN %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIAL GUESS
v0 = (zz + r.*aa).^(1-ga)/(1-ga)/rho;
v = v0;

maxit = 30;

for n=1:maxit
    V = v;
    % forward difference
    Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    Vaf(I,:) = (z + r.*amax).^(-ga); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    Vab(1,:) = (z + r.*amin).^(-ga); %state constraint boundary condition
    
    I_concave = Vab > Vaf; %indicator whether value function is concave (problems arise if this is not the case)
    
    %consumption and savings with forward difference
    cf = Vaf.^(-1/ga);
    sf = zz + r.*aa - cf;
    %consumption and savings with backward difference
    cb = Vab.^(-1/ga);
    sb = zz + r.*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = zz + r.*aa;
    Va0 = c0.^(-ga);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = sf > 0; %positive drift --> forward difference
    Ib = sb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    %make sure backward difference is used at amax
    %Ib(I,:) = 1; If(I,:) = 0;
    %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS sf > 0:
    %already taken care of automatically
    
    Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term
    
    c = Va_Upwind.^(-1/ga);
    u = c.^(1-ga)/(1-ga);
    
    %CONSTRUCT MATRIX A
    X = - min(sb,0)/da;
    Y = - max(sf,0)/da + min(sb,0)/da;
    Z = max(sf,0)/da;
    
    updiag=0; %This is needed because of the peculiarity of spdiags.
    for j=1:J
        updiag=[updiag;Z(1:I-1,j);0];
    end
    
    centdiag=reshape(Y,I*J,1);
    
    lowdiag=X(2:I,1);
    for j=2:J
        lowdiag=[lowdiag;0;X(2:I,j)];
    end
    
    AA=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);
    
    A = AA + Bswitch;
    B = (1/Delta + rho)*speye(I*J) - A;
    
    u_stacked = reshape(u,I*J,1);
    V_stacked = reshape(V,I*J,1);
    
    b = u_stacked + V_stacked/Delta;

    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = reshape(V_stacked,I,J);
    
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

ss = zz + r.*aa - c;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%
AT = A';
b = zeros(I*J,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,I*J-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg = AT\b;
g_sum = gg'*ones(I*J,1)*da*dz;
gg = gg./g_sum;

g = reshape(gg,I,J);
g_a = sum(g*dz,2);


%some graphs
figure(1)
set(gcf,'PaperPosition',[0 0 15 20])
subplot(4,2,1)
plot(a,c)
set(gca,'FontSize',14)
xlabel('Wealth, $a$','Interpreter','latex')
ylabel('$c(a,z)$','Interpreter','latex')
title('(a) Consumption, Benchmark Model','Interpreter','latex')
xlim([amin amax])
ylim([0.8 1.5])

subplot(4,2,2)
plot(a,c_robust)
set(gca,'FontSize',14)
xlabel('Wealth, $a$','Interpreter','latex')
ylabel('$c(a,z)$','Interpreter','latex')
title('(b) Consumption with Robustness Concerns','Interpreter','latex')
xlim([amin amax])
ylim([0.8 1.5])

subplot(4,2,3)
plot(a,ss,a,zeros(1,I),'k--')
set(gca,'FontSize',14)
xlabel('Wealth, $a$','Interpreter','latex')
ylabel('$s(a,z)$','Interpreter','latex')
title('(c) Saving, Benchmark Model','Interpreter','latex')
xlim([amin amax])
ylim([-0.6 0.2])

subplot(4,2,4)
plot(a,ss_robust,a,zeros(1,I),'k--')
set(gca,'FontSize',14)
xlabel('Wealth, $a$','Interpreter','latex')
ylabel('$s(a,z)$','Interpreter','latex')
title('(d) Saving with Robustness Concerns','Interpreter','latex')
xlim([amin amax])
ylim([-0.6 0.2])

subplot(4,2,5)
plot(a,g)
set(gca,'FontSize',14)
xlabel('Wealth, $a$','Interpreter','latex')
ylabel('$g(a,z)$','Interpreter','latex')
xlim([amin 2])
ylim([0 10])
title('(e) Wealth Distribution, Benchmark Model','Interpreter','latex')

subplot(4,2,6)
plot(a,g_robust)
set(gca,'FontSize',14)
xlabel('Wealth, $a$','Interpreter','latex')
ylabel('$g(a,z)$','Interpreter','latex')
xlim([amin 2])
ylim([0 10])
title('(f) Wealth Distribution with Robustness Concerns','Interpreter','latex')

subplot(4,2,7)
plot(z,h_robust)
set(gca,'FontSize',14)
xlabel('Productivity, $z$','Interpreter','latex')
ylabel('$h(a,z)$','Interpreter','latex')
xlim([zmin zmax])
title('(g) $h(a,z)$','Interpreter','latex')

subplot(4,2,8)
plot(z,drift_robust,z,mu,'k',z,zeros(J,1),'k--')
set(gca,'FontSize',14)
xlabel('Productivity, $z$','Interpreter','latex')
ylabel('$\mu(z) + h(a,z)\sigma(z)$','Interpreter','latex')
xlim([zmin zmax])
title('(h) Perceived Drift','Interpreter','latex')

print -depsc robust.eps


