% Moll lecture .. So simple that wierd result. 
% 
% Hamiltonian NK model
% GapDot = istar - r + (phi -1)*Infl
% with phi > 1
% InflDot = rho*Infl - kapa*Gap 
% SS Gap, Infl = (0,0) when istar =r
clear all
% parameter 
phi = 1.2;
dt = 0.01;
rho = 1/1.05;
ep = 1.05;
pai = 3.7;
theta =1/1.2;
kapa = (ep -1)*(1+pai)/theta;
ni=1000;
istar = ones(ni,1)*0.05 ;
r = ones(ni,1)*0.05;
% Gap = [];
% Infl = [];

%% Senario 1 
%  MIT shock in r 
% r(1) = 0.1;
% 
% Gap(1)=0;
% Infl(1) = 0;
% 
% for i= 1:ni
%   Gap(i+1) = Gap(i) + (istar(i,1) - r(i,1) + (phi -1)*Infl(i))*dt;
%   Infl(i+1) = Infl(i) + (rho*Infl(i) - kapa*Gap(i))*dt;
% end
% 
% subplot(1,2,1)
% plot(Gap)
% subplot(1,2,2)
% plot(Infl)
%%


%% Senario 2
%  permanant increase in istar
istar = istar*1.1;

Gap(1)=0;
Infl(1) = 0;
ii(1) = 0;

for i= 1:ni
  
  if istar(i,1) + phi*Infl(i) <= 0
      break
%       Gap(i+1) = Gap(i) + (-r(i,1) - Infl(i))*dt;
%       Infl(i+1) = Infl(i) + (rho*Infl(i) - kapa*Gap(i))*dt;
%       continue
  end
  Gap(i+1) = Gap(i) + (istar(i,1) - r(i,1) + (phi -1)*Infl(i))*dt;
  Infl(i+1) = Infl(i) + (rho*Infl(i) - kapa*Gap(i))*dt;
  ii(i) = istar(i,1) + phi*Infl(i);
  
end

subplot(1,3,1)
plot(Gap)
subplot(1,3,2)
plot(Infl)
subplot(1,3,3)
plot(ii)