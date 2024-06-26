%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic RBC Model 
%(Approximation in logs)
% Compute Solutions for different levels of persistence
% Introducing Government Spending Shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

%clear all;
close all;
%clear all;


%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var ly lc lk li lh ly_l lw Rk Rs Rf z g;
varexo e eg;

parameters beta chi delta alpha rho rhog sigmag gs;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
% 
% U(c,n) = log(c) + chi*log(1-n)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[param,ss] = calibgov;

alpha = param(1);
beta = param(2);
delta = param(3);
chi = param(4);
gs = ss(10); 
load parameterfile;
set_param_value('rhog',rhog)
 
rho = 0.95;
sigma   = 0.01;

sigmag = 0.01;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
  //Consumption Euler equation
  (1/exp(lc)) = beta*(1/exp(lc(+1)))*(1+alpha*(exp(lk)^(alpha-1))*exp(z(+1))*exp(lh(+1))^(1-alpha)-delta);

// Consumption/Leisure choice
 chi*exp(lc)/(1-exp(lh)) = (1-alpha)*(exp(z)*exp(lk(-1))^alpha)*(exp(lh)^(-alpha));

//Resource constraint
  exp(lc)+exp(li) + exp(g) = exp(ly);

//Production function
  exp(ly) = exp(z)*(exp(lk(-1))^alpha)*(exp(lh))^(1-alpha);

//Capital accumulation equation
  exp(li) = exp(lk)-(1-delta)*exp(lk(-1));

//Labor productivity
  exp(ly_l) = exp(ly)/exp(lh);

//Wage equation
 exp(lw) = exp(z)*(1-alpha)*exp(lk(-1))^alpha*exp(lh)^(-alpha);

//Stock return
 exp(Rs) = alpha*exp(ly)/exp(lk(-1))+1-delta;

//Capital rental rate
 exp(Rk) = alpha*exp(ly)/exp(lk(-1));

//Riskfree rate
 (1/exp(lc)) = beta*(1/exp(lc(+1)))*exp(Rf);

//TFP shock
  z = rho*z(-1)+e;

//Government spending shock
  g = (1-rhog)*log(gs) + rhog*g(-1)+eg;

end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  lk = log(ss(1));
  lc = log(ss(2));
  lh = log(ss(3));  
  li = log(ss(4));
  ly = log(ss(5));
  lw=log(ss(6));
  Rk=log(ss(7));
  Rs=log(ss(8));
  Rf=log(ss(9));
  ly_l = ly-lh;
  z = 0; 
  g = gs;
end;

shocks;
//var e = sigma^2;
var eg = sigmag^2;
end;

steady;

check;

stoch_simul(hp_filter = 1600, order = 1);
