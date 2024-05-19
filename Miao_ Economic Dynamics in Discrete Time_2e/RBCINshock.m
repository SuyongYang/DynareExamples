%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing the TFP shock with Investment Specific Technology SHock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

dynare rbc6.mod
y(:,1)=ly_e;
c(:,1)=lc_e;
n(:,1)=lh_e;
k(:,1)=lk_e;
prod(:,1)=ly_l_e;
Inv(:,1)=li_e;
Rf(:,1)=Rf_e;
Q(:,1)=Q_e;
Rs(:,1)= Rs_e;
W(:,1)=lw_e;

y(:,2)=ly_ev;
c(:,2)=lc_ev;
n(:,2)=lh_ev;
k(:,2)=lk_ev;
prod(:,2)=ly_l_ev;
Inv(:,2)=li_ev;
Rf(:,2)=Rf_ev;
Q(:,2)=Q_ev;
Rs(:,2)= Rs_ev;
W(:,2)=lw_ev;

close all;
time=1:40;
figure
subplot(3,3,1)
plot(time',y(:,1),'--',time',y(:,2),'LineWidth',2)
ylabel('Y')
legend('TFP','IST')
subplot(3,3,2)
plot(time',c(:,1),'--',time',c(:,2),'LineWidth',2)
ylabel('C')


subplot(3,3,3)
plot(time',Inv(:,1),'--',time',Inv(:,2),'LineWidth',2)
ylabel('I')

subplot(3,3,4);
plot(time',n(:,1),'--',time',n(:,2),'LineWidth',2)
ylabel('N')

subplot(3,3,5)
plot(time',prod(:,1),'--',time',prod(:,2),'LineWidth',2)
ylabel('Y/L')

subplot(3,3,6)
plot(time',W(:,1),'--',time',W(:,2),'LineWidth',2)
ylabel('W')

subplot(3,3,7);
plot(time',Q(:,1),'--',time',Q(:,2),'LineWidth',2)
ylabel('Q')

subplot(3,3,8)
plot(time',Rs(:,1),'--',time,Rs(:,2),'LineWidth',2)
ylabel('Rs')

subplot(3,3,9)
plot(time',Rf(:,1),'--',time',Rf(:,2),'LineWidth',2)
ylabel('Rf')

print -depsc figure_14_11.eps
