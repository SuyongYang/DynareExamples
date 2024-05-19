%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic RBC Model 
%(Approximation in logs)
% TFP shock and Investment-specific Technology shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

%clear all;
close all;
%clear all;


%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var ly lc lk li lh ly_l lw Rk Rs Rf z v Q;
varexo e ev;

parameters beta chi delta alpha rho rhov sigmav;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
% 
% U(c,n) = log(c) + chi*log(1-n)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[param,ss] = calibration;

alpha = param(1);
beta = param(2);
delta = param(3);
chi = param(4);
 
rho     = 0.95;  
sigma   = 0.01;

rhov=0.9;
sigmav=0.01;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
  
  //Consumption Euler equation
  exp(Q)*(1/exp(lc)) = beta*(1/exp(lc(+1)))*(exp(Rk(+1))+(1-delta)*exp(Q(+1)));

// Labor supply
 chi*exp(lc)/(1-exp(lh)) = exp(lw);

// Labor demand
 exp(lw) = exp(z)*(1-alpha)*exp(lk(-1))^alpha*exp(lh)^(-alpha);

//Resource constraint
  exp(lc)+exp(li) = exp(ly);

//Production function
  exp(ly) = exp(z)*(exp(lk(-1))^alpha)*(exp(lh))^(1-alpha);

//Capital accumulation equation
  exp(li)*exp(v) = exp(lk)-(1-delta)*exp(lk(-1));

//Capital price
  exp(Q) = 1/exp(v);

//Labor productivity
  exp(ly_l) = exp(ly)/exp(lh);

//Stock return
 exp(Rs) = (alpha*exp(ly(+1))/exp(lk)+(1-delta)*exp(Q(+1)))/exp(Q);

//Capital rental rate
 exp(Rk) = alpha*exp(ly)/exp(lk(-1));

//Riskfree rate
 (1/exp(lc)) = beta*(1/exp(lc(+1)))*exp(Rf);

//TFP shock
  z = rho*z(-1)+e;

//Investment-specific technology shock
  v = rhov*v(-1)+ev;

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
  Q = log(1);
  lw=log(ss(6));
  Rk=log(ss(7));
  Rs=log(ss(8));
  Rf=log(ss(9));
  ly_l = ly-lh;
  z = 0; 
  v = 0;
end;

shocks;
var e = sigma^2;
var ev = sigmav^2;
end;

steady;

check;

stoch_simul(order = 1);

%----------------------------------------------------------------
% 5. Some Results
%----------------------------------------------------------------

%statistic1 = 100*sqrt(diag(oo_.var(1:10,1:10)));
%dyntable('standard deviations in %',strvcat('VARIABLE','REL. S.D.'),M_.endo_names(1:10,:),statistic1,10,8,2);

%statistic2 = sqrt(diag(oo_.var(1:10,1:10)))/sqrt(diag(oo_.var(1,1)));
%dyntable('Relative standard deviations in %',strvcat('VARIABLE','REL. S.D.'),M_.endo_names(1:10,:),statistic2,10,8,2);