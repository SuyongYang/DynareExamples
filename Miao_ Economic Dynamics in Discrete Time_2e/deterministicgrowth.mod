%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deterministic Growth Model 
% Study the impact of temporary or permanent shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

%clear all;
close all;
%clear all;


%----------------------------------------------------------------
% 1. Preamble Block
%----------------------------------------------------------------

var c k n z;
varexo e;

parameters beta chi delta alpha rho ;

%----------------------------------------------------------------
% Calibration
%----------------------------------------------------------------
% 
% U(c,n) = log(c) + chi*log(1-n)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha   = 0.33;
beta    = 0.99;
delta   = 0.023;
chi     = 1.75;
rho     = 0.95;  


%----------------------------------------------------------------
% 2. Model
%----------------------------------------------------------------

model; 
  (1/c) = beta*(1/c(+1))*(1+alpha*(k^(alpha-1))*exp(z(+1))*(n(+1))^(1-alpha)-delta);
  chi*c/(1-n) = (1-alpha)*(k(-1)^alpha)*exp(z)*(n^(-alpha));
  c+ k-(1-delta)*k(-1) = (k(-1)^alpha)*exp(z)*n^(1-alpha);
  z = rho*z(-1)+e;
end;

%----------------------------------------------------------------
% 3. Steady State or Initial Values Block
%----------------------------------------------------------------

initval;
  k = 9;
  c = 0.76;
  n = 0.3;
  z = 0; 
  e = 0;
end;
steady;
check;

%---------------------------------------------
% 4. Shocks Block
%---------------------------------


%Temporary shock

%shocks;
%var e ; periods 1:10;
%values 0.01;
%end;


%Permanent shock

endval;
  k = 9;
  c = 0.76;
  n = 0.3;
  z = 0; 
  e = 0.01;
end;
steady;

shocks;
var e ; periods 1:9;
values 0;
end;


%-----------------------------------------
% 5. Computation Block
%-----------------------------------------

simul(periods=200);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Result Block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt=1:202';

subplot(2,2,1);
plot(tt,oo_.endo_simul(1,:),tt,oo_.steady_state(1)*ones(1,202),'LineWidth',2);
title('C')

subplot(2,2,2);
plot(tt,oo_.endo_simul(2,:),tt,oo_.steady_state(2)*ones(1,202),'LineWidth',2);

title('K')

subplot(2,2,3);
plot(tt,oo_.endo_simul(3,:),tt,oo_.steady_state(3)*ones(1,202),'LineWidth',2);
title('N');


subplot(2,2,4);
plot(tt,oo_.endo_simul(4,:),'LineWidth',2);
title('z')

