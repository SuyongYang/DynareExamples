clear all; close all; clc;
tic;

N=400;

Aprod_t=ones(N,1);
Aprod_t(10*5:10*10)=0.97;  %Recession hits at time t=5 for 10 years.

ga = 2; %CRRA utility with parameter gamma
d = 0.05; %depreciation rate
rho = 0.07; %discount rate
the = 0.02;
fkU = 0;
flU = 0;
fkP = 10;
flP = 8;
fyU = 0;
fyP = 0;
RS = 0.75; %important for fraction of entrepreneurs
eta = 0.4;
al = eta*RS;
be = (1-eta)*RS;
la = 3;
Aprod=1;
Bprod = 1.132;
BU = 1.3;
BP = 2.3;

%ORNSTEIN-UHLENBECK PROCESS dlog(zs) = -nu*log(z)dt + sig2*dW
%STATIONARY DISTRIBUTION IS log(z) ~ N(0,Var) WHERE Var = sig2/(2*nu)
%Corr = 0.95;
%Corr = 0.98;
%Var = 0.4;
%nu = -log(Corr);
%sig2 = 2*nu*Var;
nu = 0.02;
sig2 = 0.041/0.0512*nu;

logzmean = 0;
%zmean = exp(logzmean + Var/2); %MEAN OF LOG-NORMAL DISTRIBUTION N(0,Var)
zmin = 0.3;
zmax = 2.2;
amin = 0; %borrowing constraint
amax = 900;
%amax = 1500;

%ITERATION PARAMETERS
max_price_it = 30;
crit_price = 10^(-6);
Kcmin = 0.0000001;

rmin0 = 0;
rmax0 = 0.05;
r0 = 0.04;

J=40;
z = linspace(zmin,zmax,J);
dz = (zmax-zmin)/(J-1);
dz2 = dz^2;

I=1000;
%x=linspace(0,1,I)';
%coeff=5; power=15;
%xx=x+coeff*x.^power;
%xmax=max(xx);
%xmin=min(xx);
%a=(amax-amin)/(xmax-xmin)*xx + amin;
%daf = ones(I,1);
%dab = ones(I,1);
%daf(1:I-1)=a(2:I)-a(1:I-1);
%dab(2:I)=a(2:I)-a(1:I-1);
%daf(I)=daf(I-1); dab(1)=dab(2);

a = linspace(amin,amax,I)';
da = a(2)-a(1);
daf = ones(I,1)*da;
dab = ones(I,1)*da;

mu = (nu*(logzmean - log(z)) + sig2/2).*z; %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*z.^2; %VARIANCE (FROM ITO'S LEMMA)

aa = a*ones(1,J);
daaf = daf*ones(1,J);
daab = dab*ones(1,J);
zz = ones(I,1)*z;

maxit = 1000;
crit = 10^(-9);
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
updiag=zeros(I*J,1);
centdiag=zeros(I*J,1);
lowdiag=zeros(I*(J-1),1);

centdiag(1:I)=repmat(chi(1)+yy(1),I,1);
for j=1:J-2
    centdiag(I*j+1:I*j+I)=repmat(yy(j+1),I,1);
    lowdiag(I*(j-1)+1:j*I)=repmat(chi(j+1),I,1);
    updiag(I*j+1:I*j+I)=repmat(zeta(j),I,1);
end
centdiag((J-1)*I+1:end)=repmat(yy(J)+zeta(J),I,1);
updiag((J-1)*I+1:end)=repmat(zeta(J-1),I,1);
lowdiag(I*(J-2)+1:I*(J-1))=repmat(chi(J),I,1);

%Add up the upper, center, and lower diagonal into a sparse matrix
Bswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATION ON (w,r) STARTS HERE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmin=rmin0; rmax = rmax0; r=r0;


xi = (eta*Aprod*Bprod/(r+d)).^(1/(1-eta));
w = (1-eta)*Aprod*Bprod*xi.^eta;

%INITIAL GUESS (in iteration, use value function from previous iteration as guess)
v0 = (w.*zz + rho.*aa).^(1-ga)/(1-ga)/rho;
v = v0;

for it=1:max_price_it

r_it(it)=r;
rmin_it(it)=rmin;
rmax_it(it)=rmax;

xi = (eta*Aprod*Bprod/(r+d)).^(1/(1-eta));
w = (1-eta)*Aprod*Bprod*xi.^eta;

disp('Iteration = , r = , w = ')
disp([it,r,w])  

%FIXED COST IN UNITS OF CAPITAL
kuU = (zz.*Aprod.*BU).^(1/(1-al-be)).*(al/(r+d)).^((1-be)/(1-al-be)).*(be/w).^(be/(1-al-be)) + fkU;
kuP = (zz.*Aprod.*BP).^(1/(1-al-be)).*(al/(r+d)).^((1-be)/(1-al-be)).*(be/w).^(be/(1-al-be)) + fkP;
kU = max(min(la*aa,kuU),0);
kP = max(min(la*aa,kuP),0);
%const = la*aa<ku;
lU = (be*zz*Aprod*BU./w).^(1/(1-be)).*kU.^(al/(1-be)) + flU;
lP = (be*zz*Aprod*BP./w).^(1/(1-be)).*kP.^(al/(1-be)) + flP;
PiU = zz*Aprod*BU.*max(kU-fkU,0).^al.*max(lU-flU,0).^be - (r+d).*kU - w.*lU - fyU;
PiP = zz*Aprod*BP.*max(kP-fkP,0).^al.*max(lP-flP,0).^be - (r+d).*kP - w.*lP - fyP;
M = max(w.*zz.^the,max(PiU,PiP));
worker = (w.*zz.^the > max(PiU,PiP)); 
unproductive = (PiU > max(w.*zz.^the,PiP)); 
productive = (PiP >= max(w.*zz.^the,PiU)); 


for n=1:maxit
    V = v;
    % forward difference
    Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))./daaf(1:I-1,:);
    Vaf(I,:) = 0; %will never be used
    % backward difference
    Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))./daab(2:I,:);
    Vab(1,:) = (w*z.^the + r.*amin).^(-ga); %state constraint boundary condition
    
    I_concave = Vab > Vaf; %indicator whether value function is concave (problems arise if this is not the case)
    
    %consumption and savings with forward difference
    cf = Vaf.^(-1/ga);
    sf = M + r.*aa - cf;
    %consumption and savings with backward difference
    cb = Vab.^(-1/ga);
    sb = M + r.*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = M + r.*aa;
    Va0 = c0.^(-ga);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = sf > 0; %positive drift --> forward difference
    Ib = sb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    
    Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term
    
    c = Va_Upwind.^(-1/ga);
    u = c.^(1-ga)/(1-ga);
    
    %CONSTRUCT MATRIX A
    X = - min(sb,0)./daab;
    Y = - max(sf,0)./daaf + min(sb,0)./daab;
    Z = max(sf,0)./daaf;

    updiag=zeros(I*J,1);
    lowdiag=zeros(I*J-1,1);
    for j=1:J
	 updiag((j-1)*I+2:j*I)=Z(1:I-1,j);
         lowdiag((j-1)*I+1:j*I-1)=X(2:I,j);
    end
    centdiag=reshape(Y,I*J,1);

    B=spdiags(centdiag,0,I*J,I*J)+spdiags(updiag,1,I*J,I*J)+spdiags(lowdiag,-1,I*J,I*J);
    
    A = B + Bswitch;
    AA = (1/Delta + rho)*speye(I*J) - A;
        
    u_stacked = reshape(u,I*J,1);
    V_stacked = reshape(V,I*J,1);
    
    b = u_stacked + V_stacked/Delta;

    V_stacked = AA\b; %SOLVE SYSTEM OF EQUATIONS
    
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
da_stacked = reshape(daaf,I*J,1);
g_sum = gg'*da_stacked*dz;
gg = gg./g_sum;

g = reshape(gg,I,J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE DEMAND AND SUPPLY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
worker_s = reshape(worker,I*J,1);
Ps = reshape(productive,I*J,1);
Us = reshape(unproductive,I*J,1);
kUs = reshape(kU,I*J,1);
kPs = reshape(kP,I*J,1);
lUs = reshape(lU,I*J,1);
lPs = reshape(lP,I*J,1);
zs = reshape(zz,I*J,1);
as = reshape(aa,I*J,1);
KS(it) = (as.*gg)'*da_stacked*dz;
KUD(it) = (kUs.*Us.*gg)'*da_stacked*dz;
KPD(it) = (kPs.*Ps.*gg)'*da_stacked*dz;
KD(it) = KUD(it)+KPD(it);
LS(it) = (zs.^the.*worker_s.*gg)'*da_stacked*dz;
LUD(it) = (lUs.*Us.*gg)'*da_stacked*dz;
LPD(it) = (lPs.*Ps.*gg)'*da_stacked*dz;
LD(it) = LUD(it)+LPD(it);
Nworker(it) = (worker_s.*gg)'*da_stacked*dz;
NP(it) = (Ps.*gg)'*da_stacked*dz;
NU(it) = (Us.*gg)'*da_stacked*dz;
NE(it) = NP(it)+NU(it);

Kc = max(KS(it) - KD(it),Kcmin);

Lc(it) = Kc/xi;
ED(it) = LD(it) + Lc(it) - LS(it);

%UPDATE INTEREST RATE
if ED(it)>crit_price
    disp('ED>0, ED =, r =')
    disp([ED(it),r])
    rmax = r;
    r = 0.5*(rmax+rmin);
elseif ED(it)<-crit_price;
    disp('ED<0, ED =, r =')
    disp([ED(it),r])
    rmin = r;
    r = 0.5*(rmax+rmin);
elseif abs(ED(it))<crit_price;
    display('EQUILIBRIUM FOUND, r,w =')
    disp([r,w])
    disp('Employment Share of Corporate Sector')
    disp(Lc(it)/LS(it))
    disp('Fraction of Entrepreneurs (Productive, Unproductive) =')
    disp([NE(it),NP(it),NU(it)])
    break
end

disp('Employment Share of Corporate Sector')
disp(Lc(it)/LS(it))
disp('Fraction of Entrepreneurs (Productive, Unproductive) =')
disp([NE(it),NP(it),NU(it)])

end

v_st=V;
r_st=r;
gg_init=gg;


%----------------------------------------------------
%%%%%FIND THE TRANSITION DYNAMICS

T=100;       %The time under consideration
N=400;      %Fineness of the grid
dt=T/N;
N1=N;       %By this time, the system should have converged to new steady state

%Initial guess of the change of capital over time
r_t = r_st*ones(N,1);
rnew=r_t;       %This is just preallocation. The values will not be used.

%Preallocation
v = zeros(I,J,N);
gg_t = cell(N+1,1);
A_t=cell(N,1);
maxit = 1000;
convergence_criterion = 10^(-5);

relax=0.95;
kuU_t =zeros(I,J,N);
kuP_t =zeros(I,J,N);
kU_t = zeros(I,J,N);
kP_t = zeros(I,J,N);

lU_t = zeros(I,J,N);
lP_t = zeros(I,J,N);
PiU_t = zeros(I,J,N);
PiP_t = zeros(I,J,N);
M_t = zeros(I,J,N);
worker_t = zeros(I,J,N);
unproductive_t = zeros(I,J,N);
productive_t = zeros(I,J,N);
ED_t = zeros(N,1);

for it=1:maxit
    fprintf('ITERATION = %d\n',it);
    xi_t = (eta*Aprod_t*Bprod./(r_t+d)).^(1/(1-eta));
    w_t = (1-eta)*Aprod_t*Bprod.*xi_t.^eta;
    parfor n=1:N
        kuU_t(:,:,n) = (zz.*Aprod_t(n).*BU).^(1/(1-al-be)).*(al/(r_t(n)+d)).^((1-be)/(1-al-be)).*(be/w_t(n)).^(be/(1-al-be)) + fkU;
        kuP_t(:,:,n) = (zz.*Aprod_t(n).*BP).^(1/(1-al-be)).*(al/(r_t(n)+d)).^((1-be)/(1-al-be)).*(be/w_t(n)).^(be/(1-al-be)) + fkP;
        kU_t(:,:,n) = max(min(la*aa,kuU_t(:,:,n)),0);
        kP_t(:,:,n) = max(min(la*aa,kuP_t(:,:,n)),0);

        lU_t(:,:,n) = (be*zz*Aprod_t(n)*BU./w_t(n)).^(1/(1-be)).*kU_t(:,:,n).^(al/(1-be)) + flU;
        lP_t(:,:,n) = (be*zz*Aprod_t(n)*BP./w_t(n)).^(1/(1-be)).*kP_t(:,:,n).^(al/(1-be)) + flP;
        PiU_t(:,:,n) = zz*Aprod_t(n)*BU.*max(kU_t(:,:,n)-fkU,0).^al.*max(lU_t(:,:,n)-flU,0).^be - (r_t(n)+d).*kU_t(:,:,n) - w_t(n).*lU_t(:,:,n) - fyU;
        PiP_t(:,:,n) = zz*Aprod_t(n)*BP.*max(kP_t(:,:,n)-fkP,0).^al.*max(lP_t(:,:,n)-flP,0).^be - (r_t(n)+d).*kP_t(:,:,n) - w_t(n).*lP_t(:,:,n) - fyP;
        M_t(:,:,n) = max(w_t(n).*zz.^the,max(PiU_t(:,:,n),PiP_t(:,:,n)));
        worker_t(:,:,n) = (w_t(n).*zz.^the > max(PiU_t(:,:,n),PiP_t(:,:,n))); 
        unproductive_t(:,:,n) = (PiU_t(:,:,n) > max(w_t(n).*zz.^the,PiP_t(:,:,n))); 
        productive_t(:,:,n) = (PiP_t(:,:,n) >= max(w_t(n).*zz.^the,PiU_t(:,:,n))); 
    end
    V = v_st;
    for n=N1:-1:1
        v(:,:,n)=V;
        
        Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))./daaf(1:I-1,:);
        Vaf(I,:) = 0; %will never be used
        % backward difference
        Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))./daab(2:I,:);
        Vab(1,:) = (w_t(n)*z.^the + r_t(n).*amin).^(-ga); %state constraint boundary condition
    
        I_concave = Vab > Vaf; %indicator whether value function is concave (problems arise if this is not the case)
    
        %consumption and savings with forward difference
        cf = Vaf.^(-1/ga);
        sf = M_t(:,:,n) + r_t(n).*aa - cf;
        %consumption and savings with backward difference
        cb = Vab.^(-1/ga);
        sb = M_t(:,:,n) + r_t(n).*aa - cb;
        %consumption and derivative of value function at steady state
        c0 = M_t(:,:,n) + r_t(n).*aa;
        Va0 = c0.^(-ga);
        
        % dV_upwind makes a choice of forward or backward differences based on
        % the sign of the drift    
        If = sf > 0; %positive drift --> forward difference
        Ib = sb < 0; %negative drift --> backward difference
        I0 = (1-If-Ib); %at steady state
    
        Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term
    
        c = Va_Upwind.^(-1/ga);
        u = c.^(1-ga)/(1-ga);
    
        %CONSTRUCT MATRIX A
        X = - min(sb,0)./daab;
        Y = - max(sf,0)./daaf + min(sb,0)./daab;
        Z = max(sf,0)./daaf;
    
        updiag=zeros(I*J,1);
        lowdiag=zeros(I*J-1,1);
        for j=1:J
             updiag((j-1)*I+2:j*I)=Z(1:I-1,j);
             lowdiag((j-1)*I+1:j*I-1)=X(2:I,j);
        end
        centdiag=reshape(Y,I*J,1);
        
        A=spdiags(centdiag,0,I*J,I*J)...
            +spdiags(updiag,1,I*J,I*J)...
            +spdiags(lowdiag,-1,I*J,I*J)...
            +Bswitch;
        
        %%Note the syntax for the cell array
        A_t{n} = A;
        B = (1/dt + rho)*speye(I*J) - A;
        
        u_stacked = reshape(u,I*J,1);
        V_stacked = reshape(V,I*J,1);
        
        b = u_stacked + V_stacked/dt;
        V_stacked = B\b;
        
        V = reshape(V_stacked,I,J);
    end
    %plot(a,v(:,:,1),a,v(:,:,N))
    
    gg_t{1}=gg_init;
    for n=1:N
        AT=A_t{n}';
        %Implicit method in Updating Distribution.
        gg_t{n+1}= (speye(I*J) + AT*dt)\gg_t{n};
        %gg{n+1}=gg{n}+AT*gg{n}*dt; %This is the explicit method.
        %check(n) = gg(:,n)'*ones(2*I,1)*da;
        gg=gg_t{n};
        worker_s = reshape(worker_t(:,:,n),I*J,1);
        Ps = reshape(productive_t(:,:,n),I*J,1);
        Us = reshape(unproductive_t(:,:,n),I*J,1);
        kUs = reshape(kU_t(:,:,n),I*J,1);
        kPs = reshape(kP_t(:,:,n),I*J,1);
        lUs = reshape(lU_t(:,:,n),I*J,1);
        lPs = reshape(lP_t(:,:,n),I*J,1);
        zs = reshape(zz,I*J,1);
        as = reshape(aa,I*J,1);
        KS = (as.*gg)'*da_stacked*dz;
        KUD = (kUs.*Us.*gg)'*da_stacked*dz;
        KPD = (kPs.*Ps.*gg)'*da_stacked*dz;
        KD = KUD + KPD;
        LS = (zs.^the.*worker_s.*gg)'*da_stacked*dz;
        LUD = (lUs.*Us.*gg)'*da_stacked*dz;
        LPD = (lPs.*Ps.*gg)'*da_stacked*dz;
        LD = LUD + LPD;
        Nworker = (worker_s.*gg)'*da_stacked*dz;
        NP = (Ps.*gg)'*da_stacked*dz;
        NU = (Us.*gg)'*da_stacked*dz;
        NE = NP+NU;

        Kc = max(KS - KD,Kcmin);

        Lc = Kc/xi_t(N);
        ED_t(n) = LD + Lc - LS;
    end

    rnew=r_t;
    r_t=r_t-0.01*ED_t;

    fprintf('    Maximum change in capital is %.8f\n',max(abs(ED_t)));
    if max(abs(ED_t))<convergence_criterion
        break
    end
    
    h=figure(1);
    if mod(it,10)==0
        clf;
    end
    %plot(1:N,zeros(N,1),'r--');
    
    plot(ED_t);
    hold on;
    saveas(h,'figure.png');
    
    fig=figure(2);
    if mod(it,10)==0
        clf;
    end
    plot(r_t);
    hold on;
    saveas(fig,'interest.png');
    toc;
end
toc;
