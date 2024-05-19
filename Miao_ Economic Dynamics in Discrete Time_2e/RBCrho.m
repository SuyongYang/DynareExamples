%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing the Effects of Persistence of TFP shocks using the basic RBC
% Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

rhov=[0,0.5 0.99];
for i=1:length(rhov);
rho=rhov(i);
save parameterfile rho; 
dynare rbc2.mod noclearall
y(:,i)=ly_e;
c(:,i)=lc_e;
n(:,i)=lh_e;
k(:,i)=lk_e;
z(:,i)=z_e;
Inv(:,i)=li_e;
Rf(:,i)=Rf_e;
Rk(:,i)=Rk_e;
Rs(:,i)= Rs_e;
W(:,i)=lw_e;

end

close all;
time=1:40;
figure
subplot(2,3,1)
plot(time',z(:,1),'--',time',z(:,2),'-.', time, z(:,3),'LineWidth',2)
title('z')
legend('\rho = 0', '\rho=0.5', '\rho=0.99')

subplot(2,3,2)
plot(time',y(:,1),'--',time',y(:,2),'-.', time, y(:,3),'LineWidth',2)
title('Y')


subplot(2,3,3)
plot(time',c(:,1),'--',time',c(:,2),'-.', time, c(:,3),'LineWidth',2)
title('C')

subplot(2,3,4);
plot(time',n(:,1),'--',time',n(:,2),'-.', time, n(:,3),'LineWidth',2)
title('N')

subplot(2,3,5)
plot(time',Inv(:,1),'--',time',Inv(:,2),'-.', time, Inv(:,3),'LineWidth',2)
title('I')

subplot(2,3,6)
plot(time',k(:,1),'--',time',k(:,2),'-.', time, k(:,3),'LineWidth',2)
title('K')


print -depsc figure_14_05.eps

figure
subplot(2,2,1)
plot(time',W(:,1),'--',time',W(:,2),'-.', time, W(:,3),'LineWidth',2)
title('W')

subplot(2,2,2);
plot(time',Rk(:,1),'--',time',Rk(:,2),'-.', time, Rk(:,3),'LineWidth',2)
title('Rk')

subplot(2,2,3)
plot(time',Rs(:,1),'--',time',Rs(:,2),'-.', time, Rs(:,3),'LineWidth',2)
title('Rs')

subplot(2,2,4)
plot(time',Rf(:,1),'--',time',Rf(:,2),'-.', time, Rf(:,3),'LineWidth',2)
title('Rf')
legend('\rho = 0', '\rho=0.5', '\rho=0.99')
print -depsc figure_14_06.eps
