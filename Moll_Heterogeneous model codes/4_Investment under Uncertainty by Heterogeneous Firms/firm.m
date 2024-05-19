%Thanks to Peter Tian for suggesting this problem
clear all; close all; clc;
tic;

rho = 0.05; %discount rate  %% attention

%ORNSTEIN-UHLENBECK PROCESS dlog(z) = -the*log(z)dt + sig2*dW %%%  here
% sig2 should be sig,which represents standard deviation
%STATIONARY DISTRIBUTION IS log(z) ~ N(0,Var) WHERE Var = sig2/(2*the)
Var = 0.026^2;  %%% z corresponds to the  idiosyncratic shock in equation 3
zmean = exp(Var/2); 
Corr = 0.859;
the = -log(Corr);
sig2 = 2*the*Var;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=40;
zmin = zmean*0.6;
zmax = zmean*1.4;
z = linspace(zmin,zmax,J);
dz = (zmax-zmin)/(J-1);
dz2 = dz^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=100;
kmin = 1;
kmax = 100;
k = linspace(kmin,kmax,I)';
dk = (kmax-kmin)/(I-1);
%plot(z,lognpdf(z,0,Var)) %PLOT STATIONARY DISTRIBUTION AND CHECK THAT CHOICE OF GRID DOESN'T CUT OFF TOO MUCH OF TAILS

kk = k*ones(1,J);
zz = ones(I,1)*z;

mu = (-the*log(z) + sig2/2).*z; 
s2 = sig2*(z.^2);
maxit= 20;
crit = 10^(-6);
Delta = 1000;
%%%%%%%%%%set parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha= 0.5;     %curvature in production function
theta=2.7;%%%% quadratic adjustment cost
delta=0.025; %%%depreciation rate
F = zz.*kk.^alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vkf = zeros(I,J);
Vkb = zeros(I,J);
Vzf = zeros(I,J);
Vzb = zeros(I,J);
Vzz = zeros(I,J);
c = zeros(I,J);  %%c represents xestment in equation 3

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

%INITIAL GUESS
v0 =(F - delta*kk - 0.5*theta*delta^2*kk)/rho;    %%%% attention please, need carefully check. the U value at  steady state 
v = v0;

plot(k,v)

for n=1:maxit
    V = v;
    % forward difference
    Vkf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/dk;
    Vkf(I,:) =   (1+theta*delta)*ones(1,J); %state constraint boundary condition
    % backward difference
    Vkb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/dk;
    Vkb(1,:) =   (1+theta*delta)*ones(1,J); %state constraint boundary condition
    
    %xestment with forward difference
    xf = (Vkf-1)/theta.*kk;
    sf = xf-delta.*kk;
    Hf = F - xf - 0.5*theta*(xf./kk).^2.*kk + Vkf.*sf;
    
    %xestment with backward difference
    xb = (Vkb-1)/theta.*kk;
    sb =  xb-delta.*kk; 
    Hb = F - xb - 0.5*theta*(xb./kk).^2.*kk + Vkb.*sb;
    
    %xestment at steady state
    x0 =delta*kk;
    
    Ineither = (1-(sf>0)) .* (1-(sb<0));
    Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
    Iboth = (sb<0).*(sf>0);
    Ib = Iunique.*(sb<0) + Iboth.*(Hb>=Hf);
    If = Iunique.*(sf>0) + Iboth.*(Hf>=Hb);
    I0 = Ineither;

    x = xf.*If + xb.*Ib + x0.*I0;
    profits = F - x - 0.5*theta*(x./kk).^2.*kk;
    
    %CONSTRUCT MATRIX A
    X = - sb.*Ib/dk;
    Y = - sf.*If/dk + sb.*Ib/dk;
    Z = sf.*If/dk;
    
    %%%%spdiags   when updiag,0 first;when lowdiag,0 last.
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
    
    profits_stacked = reshape(profits,I*J,1);
    V_stacked = reshape(V,I*J,1);
    
    b = profits_stacked + V_stacked/Delta;

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

plot(k,v)

plot(dist)

kdot = x - delta.*kk;
plot(k,kdot,k,zeros(I,1))


