%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing the Effects of Different Values of Risk Aversion Parameter Using the Basic RBC
% Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

gammav=[0.75, 1, 2, 10];
for i=1:length(gammav);
gamma = gammav(i)
save parameterfile gamma; 
dynare rbc3.mod noclearall
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
subplot(3,3,1)
plot(time',y(:,1),'--',time',y(:,2),':', time, y(:,3),time, y(:,4),'-.','LineWidth',2)
title('Y')


subplot(3,3,2)
plot(time',c(:,1),'--',time', c(:,2),':', time, c(:,3),time, c(:,4),'-.','LineWidth',2)
title('C')

subplot(3,3,3);
plot(time',n(:,1),'--',time',n(:,2),':', time, n(:,3),time,n(:,4),'-.','LineWidth',2)
title('N')

subplot(3,3,4)
plot(time',Inv(:,1),'--',time',Inv(:,2),':', time, Inv(:,3),time,Inv(:,4),'-.', 'LineWidth',2)
title('I')

subplot(3,3,5)
plot(time',k(:,1),'--',time',k(:,2),':', time, k(:,3),time, k(:,4),'-.', 'LineWidth',2)
title('K')

subplot(3,3,6)
plot(time',W(:,1),'--',time',W(:,2),':', time, W(:,3),time,W(:,4),'-.', 'LineWidth',2)
title('W')

subplot(3,3,7);
plot(time',Rk(:,1),'--',time',Rk(:,2),':', time, Rk(:,3),time, Rk(:,4), '-.', 'LineWidth',2)
title('Rk')

subplot(3,3,8)
plot(time',Rs(:,1),'--',time',Rs(:,2),':', time, Rs(:,3),time, Rs(:,4), '-.', 'LineWidth',2)
title('Rs')

subplot(3,3,9)
plot(time',Rf(:,1),'--',time',Rf(:,2),':', time, Rf(:,3), time, Rf(:,4), '-.', 'LineWidth',2)
title('Rf')
legend('\gamma=0.75', '\gamma = 1', '\gamma=5', '\gamma=10')
print -depsc figure_14_07.eps

