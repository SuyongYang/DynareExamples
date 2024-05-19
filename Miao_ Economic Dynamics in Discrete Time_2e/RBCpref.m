%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing the Effects of Persistence of Preference shocks using the basic RBC
% Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
phiv=[0,0.5 0.95];
for i=1:length(phiv);
rhophi = phiv(i);
save parameterfile rhophi; 
dynare rbc8.mod noclearall
y(:,i)=ly_ephi;
c(:,i)=lc_ephi;
n(:,i)=lh_ephi;
k(:,i)=lk_ephi;
phi(:,i)=phi_ephi;
Inv(:,i)=li_ephi;
Rf(:,i)=Rf_ephi;
Rk(:,i)=Rk_ephi;
Rs(:,i)= Rs_ephi;
W(:,i)=lw_ephi;

end

close all;
time=1:40;
figure
subplot(3,3,1)
plot(time',phi(:,1),'--',time',phi(:,2),'-.', time, phi(:,3),'LineWidth',2)
title('Preference Shock')
legend('\rho_p = 0', '\rho_p=0.5', '\rho_p=0.95')

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
plot(time',Rs(:,1),'--',time',Rs(:,2),'-.', time, Rk(:,3),'LineWidth',2)
title('Rs')


subplot(3,3,9)
plot(time',Rf(:,1),'--',time',Rf(:,2),'-.', time, Rf(:,3),'LineWidth',2)
title('Rf')


print -depsc figure_14_13.eps
