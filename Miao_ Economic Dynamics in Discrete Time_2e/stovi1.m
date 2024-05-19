clear all
%
% Solve the stochastic OGM by value iteration (brute force)
%
gamma	= 1.5;
delta	= 0.1;
beta	= 0.95;
alpha	= 0.30;
p		= 0.9;
PI		= [p 1-p;1-p p];
se		= 0.2;
ab		= 0;
am		= exp(ab-se);
as		= exp(ab+se);
A=[am as];
nba=2;

ks		= ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
csy	= 1-alpha*beta*delta/(1-beta*(1-delta));
dev	= 0.9;
kmin	= (1-dev)*ks;
kmax	= (1+dev)*ks;
nbk	= 100;
devk	= (kmax-kmin)/(nbk-1);
k		= linspace(kmin,kmax,nbk)';
c		= zeros(nbk,nba);
u		= zeros(nbk,nba);
v		= zeros(nbk,nba);
Tv		= zeros(nbk,nba);
iter	= 1;
crit	= 1;
tol	= 1e-6;
while crit>tol;
   for i=1:nbk
      for j=1:nba;
         c	= A(j)*k(i)^alpha+(1-delta)*k(i)-k;
         neg		= find(c<=0);     
         c(neg)	= NaN;
         u(:,j)	= (c.^(1-gamma)-1)/(1-gamma);
         u(neg,j)	= -1e12;
      end
      [Tv(i,:),dr(i,:)]	= max(u+beta*(v*PI));
   end;
   crit	= max(max(abs(Tv-v)));
   v		= Tv;
   iter	= iter+1;
end

kp		= k(dr);
for j=1:nba;
   c(:,j)	= A(j)*k.^alpha+(1-delta)*k-kp(:,j);
end
u		= (c.^(1-gamma)-1)/(1-gamma);
v		= u/(1-beta);

figure

subplot(121);
h=plot(k,[k kp]);
title('Capital policy function')
xlabel('k')
subplot(122);
h=plot(k,c);
title('Consumption policy function')
xlabel('k')
figure
plot(k,v);
title('Value function')
xlabel('k')
%print -depsc2 valuevi.eps
