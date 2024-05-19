figure

subplot(221)
icut = 60;
acut =  results.a(1:icut);
vcut = results.V(1:icut,:);
z    = results.z;

set(gca,'FontSize',14)
mesh(acut,z,vcut')
view([45 25])
xlabel('Wealth, $a$','FontSize',14,'interpreter','latex')
ylabel('Productivity, $z$','FontSize',14,'interpreter','latex')
xlim([parameters.amin max(acut)])
ylim([parameters.zmin parameters.zmax])
title('Value function $j(a,z)$','FontSize',14,'interpreter','latex')

subplot(222)
icut = 60;
acut = results.a(1:icut);
ccut = results.c(1:icut,:);

set(gca,'FontSize',14)
mesh(acut,z,ccut')
view([45 25])
xlabel('Wealth, $a$','FontSize',14,'interpreter','latex')
ylabel('Productivity, $z$','FontSize',14,'interpreter','latex')
xlim([parameters.amin max(acut)])
ylim([parameters.zmin parameters.zmax])
title('Consumption $c(a,z)$','FontSize',14,'interpreter','latex')

subplot(223)
icut = 60;
acut = results.a(1:icut);
sscut = results.s(1:icut,:);

set(gca,'FontSize',14)
mesh(acut,z,sscut')
view([45 25])
xlabel('Wealth, $a$','FontSize',14,'interpreter','latex')
ylabel('Productivity, $z$','FontSize',14,'interpreter','latex')
xlim([parameters.amin max(acut)])
ylim([parameters.zmin parameters.zmax])
title('Savings $s(a,z)$','FontSize',14,'interpreter','latex')

subplot(224)
gcut = results.g(1:icut,:);
set(gca,'FontSize',14)
mesh(acut,z,gcut')
view([45 25])
xlabel('Wealth, $a$','FontSize',14,'interpreter','latex')
ylabel('Productivity, $z$','FontSize',14,'interpreter','latex')
xlim([parameters.amin max(acut)])
ylim([parameters.zmin parameters.zmax])
zlim([-0.01 0.7])
title('Wealth-productivity density $g(a,z)$','FontSize',14,'interpreter','latex')

%%
% figure
% subplot(121)
% icut = 60;
% acut = results.a(1:icut);
% ccut = results.c(1:icut,:);
% 
% set(gca,'FontSize',14)
% mesh(acut,z,ccut')
% view([45 25])
% xlabel('Wealth, $a$','FontSize',14,'interpreter','latex')
% ylabel('Productivity, $z$','FontSize',14,'interpreter','latex')
% xlim([parameters.amin max(acut)])
% ylim([parameters.zmin parameters.zmax])
% title('Consumption $c(a,z)$','FontSize',14,'interpreter','latex')
% 
% 
% 
% subplot(122)
% gcut = results.g(1:icut,:);
% set(gca,'FontSize',14)
% mesh(acut,z,gcut')
% view([45 25])
% xlabel('Wealth, $a$','FontSize',14,'interpreter','latex')
% ylabel('Productivity, $z$','FontSize',14,'interpreter','latex')
% xlim([parameters.amin max(acut)])
% ylim([parameters.zmin parameters.zmax])
% zlim([-0.01 0.7])
% title('Wealth-productivity distribution $g(a,z)$','FontSize',14,'interpreter','latex')
