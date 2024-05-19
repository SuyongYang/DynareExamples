%Optimized for speed by SeHyoun Ahn

clear all; close all; clc;
tic;
ga = 2; %CRRA utility with parameter gamma
rho = 0.05; %discount rate
alpha = 0.3; %CURVATURE OF PRODUCTION FUNCTION
d = 0.05; %DEPRECIATION RATE


%ORNSTEIN-UHLENBECK PROCESS dlog(z) = -the*log(z)dt + sig*dW
%STATIONARY DISTRIBUTION IS log(z) ~ N(0,Var) WHERE Var = sig2/(2*the)
Var = 0.07;
zmean = exp(Var/2); %MEAN OF LOG-NORMAL DISTRIBUTION N(0,Var)
Corr = 0.9;
the = -log(Corr);
sig2 = 2*the*Var;

k_st = (alpha*zmean/(rho+d))^(1/(1-alpha));

I=100;
kmin = 0.3*k_st;
kmax = 3*k_st;
k = linspace(kmin,kmax,I)';
dk = (kmax-kmin)/(I-1);

J=40;
zmin = zmean*0.8;
zmax = zmean*1.2;
z = linspace(zmin,zmax,J);
dz = (zmax-zmin)/(J-1);
dz2 = dz^2;

%plot(z,lognpdf(z,0,Var)) %PLOT STATIONARY 
%DISTRIBUTION AND CHECK THAT CHOICE OF GRID 
%DOESN'T CUT OFF TOO MUCH OF TAILS

kk = k*ones(1,J);
zz = ones(I,1)*z;

%behavior of z : dz = mu*dt + s*dW
mu = (-the*log(z) + sig2/2*z.^(-2)).*z;   %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*z.^2;   %VARIANCE (FROM ITO'S LEMMA)

maxit= 100;
crit = 10^(-6);
Delta = 1000;

Vaf = zeros(I,J);
Vab = zeros(I,J);
Vzf = zeros(I,J);
Vzb = zeros(I,J);
Vzz = zeros(I,J);
c = zeros(I,J);

%CONSTRUCT MATRIX Bswitch SUMMARIZING EVOLUTION OF z
yy = - s2/dz2 - mu/dz;
chi =  s2/(2*dz2);
zeta = mu/dz + s2/(2*dz2);

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


%INITIAL GUESS
v0 = (zz.*kk.^alpha).^(1-ga)/(1-ga)/rho;
v = v0;

maxit = 30;

for n=1:maxit
    V = v;
    % forward difference
    Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/dk;
    Vaf(I,:) = (z.*kmax.^alpha - d.*kmax).^(-ga); %will never be used, but impose state constraint a<=kmax just in case
    % backward difference
    Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))/dk;
    Vab(1,:) = (z.*kmin.^alpha - d.*kmin).^(-ga); %will never be used, but impose state constraint a<=kmax just in case
    
    I_concave = Vab > Vaf; %indicator whether value function is concave (problems arise if this is not the case)
    
    %consumption and savings with forward difference
    cf = Vaf.^(-1/ga);
    sf = zz.*kk.^alpha - d.*kk - cf;
    %consumption and savings with backward difference
    cb = Vab.^(-1/ga);
    sb = zz.*kk.^alpha - d.*kk - cb;
    %consumption and derivative of value function at steady state
    c0 = zz.*kk.^alpha - d.*kk;
    Va0 = c0.^(-ga);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = sf > 0; %positive drift --> forward difference
    Ib = sb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    %make sure backward difference is used at kmax
    %Ib(I,:) = 1; If(I,:) = 0;
    %STATE CONSTRAINT at kmin: USE BOUNdkRY CONDITION UNLESS sf > 0:
    %already taken care of automatically
    
    Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; 
    %important to include third term
    
    c = Va_Upwind.^(-1/ga);
    u = c.^(1-ga)/(1-ga);
    
    %CONSTRUCT MATRIX A
    X = - min(sb,0)/dk;
    Y = - max(sf,0)/dk + min(sb,0)/dk;
    Z = max(sf,0)/dk;
    
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
    
    b = u_stacked + V_stacked/Delta; %implicit method

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

ss = zz.*kk.^alpha - d.*kk - c;
set(gca,'FontSize',14)
plot(k,ss,k,zeros(1,I),'--')
xlabel('k')
ylabel('s(k,z)')
xlim([kmin kmax])