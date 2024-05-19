% Use Farmer, Khramov, and Nicolo (2015) Method to solve a DNK model
% The determinate case
% 
%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;


%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var x R pi zx zpi;
varexo ex epi eR;

parameters beta gamma kappa rhoR rhox rhopi psi1 psi2;


beta    = 0.99;
gamma   = 1;
kappa   = 0.77;  
rhoR    = 0.95;
rhox    = 0.68;
rhopi   = 0.82;
psi1    = 1.1;
psi2    = 0.17; 

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model(linear); 
//(1) IS curve
 x = x(+1)-1/gamma*(R-pi(+1))+zx;
//(2) NKPC
 pi = beta*pi(+1)+kappa*x+zpi;
//(3) Interest rate rule
 R = rhoR*R(-1)+(1-rhoR)*(psi1*pi+psi2*x)+eR;
//(4) Shock to IS curve
  zx = rhox*zx(-1)+ex;
//(5) Shock to NKPC
 zpi = rhopi*zpi(-1)+epi;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  x = 0;
  pi = 0;
  R = 0;
  zx = 0; 
  zpi = 0;
end;

shocks;
var eR = 0.01/4;
var ex = 0.01;
var epi = 0.01/4;
end;

check;

stoch_simul;
