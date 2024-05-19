%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing the Effects of Persistence of Government Spending shocks using the basic RBC
% Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
rhogv=[0,0.5 0.95];
for i=1:length(rhogv);
rhog = rhogv(i);
save parameterfile rhog; 
dynare rbc7.mod noclearall
y(:,i)=ly_eg;
c(:,i)=lc_eg;
n(:,i)=lh_eg;
k(:,i)=lk_eg;
g(:,i)=g_eg;
Inv(:,i)=li_eg;
Rf(:,i)=Rf_eg;
Rk(:,i)=Rk_eg;
Rs(:,i)= Rs_eg;
W(:,i)=lw_eg;

end

close all;
time=1:40;
figure
subplot(3,3,1)
plot(time',g(:,1),'--',time',g(:,2),'-.', time, g(:,3),'LineWidth',2)
title('G')
legend('\rho_g = 0', '\rho_g=0.5', '\rho_g=0.95')

subplot(3,3,2)
plot(time',y(:,1),'--',time',y(:,2),'-.', time, y(:,3),'LineWidth',2)
title('Y')


subplot(3,3,3)
plot(time',c(:,1),'--',time',c(:,2),'-.', time, c(:,3),'LineWidth',2)
title('C')

subplot(3,3,4);
plot(time',n(:,1),'--',time',n(:,2),'-.', time, n(:,3),'LineWidth',2)
title('N')

subplot(3,3,5)
plot(time',Inv(:,1),'--',time',Inv(:,2),'-.', time, Inv(:,3),'LineWidth',2)
title('I')

subplot(3,3,6)
plot(time',k(:,1),'--',time',k(:,2),'-.', time, k(:,3),'LineWidth',2)
title('K')

%figure
subplot(3,3,7)
plot(time',W(:,1),'--',time',W(:,2),'-.', time, W(:,3),'LineWidth',2)
title('W')

subplot(3,3,8);
plot(time',Rk(:,1),'--',time',Rk(:,2),'-.', time, Rk(:,3),'LineWidth',2)
title('Rk')


subplot(3,3,9)
plot(time',Rf(:,1),'--',time',Rf(:,2),'-.', time, Rf(:,3),'LineWidth',2)
title('Rf')


print -depsc figure_14_12.eps
