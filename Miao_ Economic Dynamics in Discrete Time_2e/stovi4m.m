
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solve the stochastic OGM by value function iteration method
%
%    V(K,A) = max (C^(1-gamma)-1)/(1-gamma)+ beta*E[V(K',A')]
%      subject to C+K'=A*K^alpha+(1-delta)*K
%      A is a Markov process
%
% Use interpolation of the value function
% Written by Jianjun Miao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

clear all;
close all

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma	= 1.5;
delta	= 0.1;
beta	= 0.95;
alpha	= 0.30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p		= 0.9;
PI		= [p 1-p;1-p p];
se		= 0.2;
ab		= 0;
am		= exp(ab-se);
as		= exp(ab+se);
A=[am as];

nba=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deterministic Steady State
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ks		= ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
csy	= 1-alpha*beta*delta/(1-beta*(1-delta));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbk	= 50;   % less grid points for the state variable (capital)
nbc	= 2000; % more grid points for the choice variable, which is outside the capital grid
dev	= 0.9;

kmin	= ks/2; %
kmax	=ks*2;
k		= linspace(kmin,kmax,nbk)';  %Capital grid
kpol	= linspace(kmin,kmax,nbc)';  %Put more grid points for Choice of the next period capital 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Period utility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

util	= zeros(nbc,nbk,nba);

for j=1:nba;
   cons		  	= repmat(A(j)*k'.^alpha+(1-delta)*k',nbc,1)-repmat(kpol,1,nbk);      
   util(:,:,j) = (cons.^(1-gamma)-1)/(1-gamma);
   utilj 		= util(:,:,j);
   neg		   = find(cons<=0);  % ensure power function is well defined 
   utilj(neg)  = -inf;           % never choose negative consumption
   util(:,:,j) = utilj;
  clear cons neg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Core Value Function Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v		= zeros(nbk,nba);
Tv		= zeros(nbk,nba);
iter	= 1;
crit	= 1;
tol	= 1e-6;

while crit>tol;

      for j=1:nba;       
         vi		       = interp1(k,v,kpol,'spline'); %nbc by nba matrix
         [vtmp,drtmp]  = max(util(:,:,j)+beta*repmat(vi*PI(j,:)',1,nbk));
         Tv(:,j)      = vtmp';
         dr(:,j)      = drtmp';
      end
  
   crit	= max(max(abs(Tv-v)));
   v		= Tv;
   disp(sprintf('Iteration # %2d \tCriterion: %g',iter,crit))
   iter	= iter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consumption and Capital Policy functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kp=zeros(nbk,nba);
ct=zeros(nbk,nba);
for j=1:nba;
 kp(:,j)=kpol(dr(:,j));
 ct(:,j)	= A(j)*k.^alpha+(1-delta)*k-kp(:,j);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%
% Figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
figure

subplot(131);

plot(k,kp(:,1),k,kp(:,2),'--',k,k,'-.','LineWidth',1.5);
xlabel('$K$','interpreter','latex','fontsize',12)
title('$K^{\prime}=g(K,A)$','interpreter','latex','fontsize',12)

legend({'$A=A_l$','$A=A_h$','$45^{0}$-line'},'interpreter','latex','location','southeast','fontsize',12)

subplot(132);plot(k,ct(:,1),k,ct(:,2),'--','linewidth',1);
xlabel('$K$','interpreter','latex','fontsize',12)
title('$C(K,A)$','interpreter','latex','fontsize',12)

subplot(133)
plot(k,v(:,1),k,v(:,2),'--','linewidth',1.5);
xlabel('$K$','interpreter','latex','fontsize',12)
title('$V(K,A)$','interpreter','latex','fontsize',12)

%pause
%print -depsc2 vfstovi.eps
 
%close