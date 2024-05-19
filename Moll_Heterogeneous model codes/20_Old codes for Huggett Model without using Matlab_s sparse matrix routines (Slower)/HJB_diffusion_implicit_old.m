clear all; close all; clc;
tic;
ga = 2; %CRRA utility with parameter gamma
r = 0.03; %interest rate
rho = 0.05; %discount rate

%ORNSTEIN-UHLENBECK PROCESS dlog(z) = -the*log(z)dt + sig2*dW
%STATIONARY DISTRIBUTION IS log(z) ~ N(0,Var) WHERE Var = sig2/(2*the)
Var = 0.07;
zmean = exp(Var/2); %MEAN OF LOG-NORMAL DISTRIBUTION N(0,Var)
Corr = 0.9;
the = -log(Corr);
sig2 = 2*the*Var;


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

mu = (-the*log(z) + sig2/2).*z; %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*z.^2; %VARIANCE (FROM ITO'S LEMMA)

maxit= 100;
crit = 10^(-6);
Delta = 1000;

Vaf = zeros(I,J);
Vab = zeros(I,J);
Vzf = zeros(I,J);
Vzb = zeros(I,J);
Vzz = zeros(I,J);
c = zeros(I,J);

%CONSTRUCT MATRIX Aswitch SUMMARIZING EVOLUTION OF z
yy = - s2/dz2 - mu/dz;
chi =  s2/(2*dz2);
zeta = mu/dz + s2/(2*dz2);

Aswitch = zeros(I*J,I*J);
%CASE 2<=j<=J-1
for j=2:J-1
    for i=1:I
        row = (j-1)*I+i;
        i_j = (j-1)*I+i; %index for (i,j) entry
        i_jp1 = j*I+i; %index for (i,j+1) entry
        i_jm1 = (j-2)*I+i; %index for (i,j-1) entry
        Aswitch(row,i_j) = yy(j);
        Aswitch(row,i_jm1) = chi(j);
        Aswitch(row,i_jp1) = zeta(j);
    end
end

%CASE j=1
j=1;
for i=1:I
    row = (j-1)*I+i;
    i_j = (j-1)*I+i; %index for (i,j) entry
    i_jp1 = j*I+i; %index for (i,j+1) entry
    Aswitch(row,i_j) = yy(j) + chi(j);
    Aswitch(row,i_jp1) = zeta(j);
end
   
%CASE j=J
j=J;
for i=1:I
    row = (j-1)*I+i;
    i_j = (j-1)*I+i; %index for (i,j) entry
    i_jp1 = j*I+i; %index for (i,j+1) entry
    i_jm1 = (j-2)*I+i; %index for (i,j-1) entry
    Aswitch(row,i_j) = yy(j) + zeta(j);
    Aswitch(row,i_jm1) = chi(j);
end

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
    
    for j=1:J
        for i=2:I-1
            AAj(i,i-1) = X(i,j);
            AAj(i,i) = Y(i,j);
            AAj(i,i+1) = Z(i,j);
        end
        AAj(1,1)=Y(1,j); AAj(1,2) = Z(1,j);
        AAj(I,I)=Y(I,j); AAj(I,I-1) = X(I,j);
        %create block-diagonal matrix
        AA((j-1)*I+1:j*I,(j-1)*I+1:j*I) = AAj;
    end
   
    A = AA + Aswitch;
    B = (rho + 1/Delta)*eye(I*J) - A;
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
set(gca,'FontSize',14)
plot(a,ss,a,zeros(1,I),'--')
xlabel('a')
ylabel('s(a,z)')
xlim([amin amax])