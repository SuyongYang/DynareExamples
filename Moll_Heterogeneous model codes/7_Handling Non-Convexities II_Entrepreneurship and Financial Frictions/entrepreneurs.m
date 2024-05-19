clear all; close all; clc;
tic;
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
% plot(z,lognpdf(z,logzmean,sqrt(Var)))
% ylim([0 1.5])

I=1000;
% x = linspace(0,1,I)';
% coeff = 5; power = 15;
% xx  = x + coeff*x.^power;
% xmax = max(xx); xmin = min(xx);
% a = (amax-amin)/(xmax - xmin)*xx + amin;
% daf = ones(I,1);
% dab = ones(I,1);
% daf(1:I-1) = a(2:I)-a(1:I-1);
% dab(2:I) = a(2:I)-a(1:I-1);
% daf(I)=daf(I-1); dab(1)=dab(2);

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

maxit = 20;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATION ON (w,r) STARTS HERE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmin=rmin0; rmax = rmax0; r=r0;

%max_price_it = 1;
for it=1:max_price_it

r_it(it)=r;
rmin_it(it)=rmin;
rmax_it(it)=rmax;

xi = (eta*Aprod*Bprod/(r+d)).^(1/(1-eta));
w = (1-eta)*Aprod*Bprod*xi.^eta;

disp('Iteration = , r = , w = ')
disp([it,r,w])

%INITIAL GUESS (in iteration, use value function from previous iteration as guess)
v0 = (w.*zz + rho.*aa).^(1-ga)/(1-ga)/rho;
v = v0;
  

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
    
    updiag=0; %This is needed because of the peculiarity of spdiags.
    for j=1:J
        updiag=[updiag;Z(1:I-1,j);0];
    end
    
    centdiag=reshape(Y,I*J,1);
    
    lowdiag=X(2:I,1);
    for j=2:J
        lowdiag=[lowdiag;0;X(2:I,j)];
    end
    
    B=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);
    
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

%WORK WITH RESCALED DENSITY \tilde{g} BELOW
da_tilde = 0.5*(dab + daf);
da_tilde(1) = 0.5*daf(1); da_tilde(I) = 0.5*dab(I);
da_stacked = reshape(da_tilde*ones(1,J),I*J,1);
grid_diag = spdiags(da_stacked,0,I*J,I*J);

AT = A';
b = zeros(I*J,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,I*J-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
g_tilde = AT\b;
%rescale \tilde{g} so that it sums to 1
g_sum = g_tilde'*ones(I*J,1);
g_tilde = g_tilde./g_sum;

gg = (dz*grid_diag)\g_tilde; %convert from \tilde{g} to g

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



plot([rmin_it;rmax_it;r_it]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE EQUILIBRIUM QUANTITIES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SAVING POLICY FUNCTION
ss = M + r.*aa - c;

%WEALTH DISTRIBUTION
g_a = sum(g.*dz,2);
%CDF
for i=1:I
    G_a(i)=sum(g_a(1:i).*daf(1:i));
end

p = linspace(0,1,100);
for pp=1:100
    [val index(pp)] = min(abs(G_a - p(pp)));
end
perc = a(index);

%WEALTH SHARES
integrand = a.*g_a.*daf;
total = sum(integrand);
bottom50 = sum(integrand(1:index(50)-1))/total*100;
next40 = sum(integrand(index(50):index(90)-1))/total*100;
next9 = sum(integrand(index(90):index(99)-1))/total*100;
top1 = sum(integrand(index(99):I))/total*100;

check = bottom50 + next40 + next9 + top1
disp('bottom50,next40,next9,top1 =')
disp([bottom50,next40,next9,top1])

disp('50,90,99 percentile')
disp([a(index(50)),a(index(90)),a(index(99))])


%SHARE OF WEALTH HELD BY ENTREPRENEURS
prod_wealth_share = sum((as.*Ps.*gg)'*da_stacked*dz)/total*100;
unprod_wealth_share = sum((as.*Us.*gg)'*da_stacked*dz)/total*100;
ent_wealth_share = prod_wealth_share+unprod_wealth_share;
worker_wealth_share = sum((as.*worker_s.*gg)'*da_stacked*dz)/total*100;
ent_wealth_share + worker_wealth_share;
disp('Fraction of Entrepreneurs (Productive, Unproductive) =')
disp([NE(it),NP(it),NU(it)])
disp('Share of Wealth held by Entrepreneurs (Productive, Unproductive) =')
disp([ent_wealth_share,prod_wealth_share,unprod_wealth_share])


%%%%%%%%%%
% GRAPHS %
%%%%%%%%%%

jmed = 32; jhigh = 35;
amax_fig = 80;

%CONSUMPTION AND SAVING POLICY FUNCTIONS
close all
figure(1)
set(gcf,'PaperPosition',[0 0 15 5])
subplot(1,2,1)
set(gca,'FontSize',16)
plot(a,[ss(:,1),ss(:,jmed),ss(:,jhigh)],a,zeros(1,I),'k--','LineWidth',2)
legend('low z','medium z','high z','Location','NorthWest')
xlabel('Wealth, a','FontSize',16,'interpreter','latex')
ylabel('Saving $s(a,z)$','FontSize',16,'interpreter','latex')
title('(a) Saving, with Entrepreneurship','interpreter','latex')
xlim([amin amax_fig])
ylim([-2 7])

subplot(1,2,2)
set(gca,'FontSize',16)
plot(a,[c(:,1),c(:,jmed),c(:,jhigh)],'LineWidth',2)
legend('low z','medium z','high z','Location','NorthWest')
xlabel('Wealth, a','FontSize',16,'interpreter','latex')
ylabel('Consumption $c(a,z)$','FontSize',16,'interpreter','latex')
title('(c) Consumption, with Entrepreneurship','interpreter','latex')
xlim([amin amax_fig])
ylim([0.5 8])
print -depsc fig_policies.eps

%WEALTH DISTRIBUTION
figure(2)
set(gca,'FontSize',16)
plot(a,g_a,'LineWidth',2)
xlabel('Wealth, a','FontSize',16,'interpreter','latex')
ylabel('Density','FontSize',16,'interpreter','latex')
title('(a) with Entrepreneurship','interpreter','latex')
xlim([amin amax_fig])
ylim([0 0.06])
print -depsc fig_wealthdist.eps
