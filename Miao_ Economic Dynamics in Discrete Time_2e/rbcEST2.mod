%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic RBC Model 
%(Approximation in logs)
% Generate Simulated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------
% 
% U(c,n) = log(c) + chi*log(1-n)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------


close all;
clc;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var ly lc lk li lh lw Rk z;
varexo e;

parameters beta delta chi alpha rho;

%beta = 0.98;
%delta = 0.025;
%chi =1.75;
%alpha = 0.3;
%rho =0.9;


%----------------------------------------------------------------
% 2. Declaring observable variables
%---------------------------------------------------------------

varobs ly;


%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
  
  //1. Consumption Euler equation
 1/exp(lc) = beta*(1/exp(lc(+1)))*(1+alpha*(exp(lk)^(alpha-1))*exp(z(+1))*exp(lh(+1))^(1-alpha)-delta);

//2. Labor supply
 chi*exp(lc)/(1-exp(lh)) = exp(lw);

//3. Labor demand
 exp(lw) = exp(z)*(1-alpha)*exp(lk(-1))^alpha*exp(lh)^(-alpha);

//3. Resource constraint
  exp(lc)+exp(li) = exp(ly);

//5. Production function
  exp(ly) = exp(z)*(exp(lk(-1))^alpha)*(exp(lh))^(1-alpha);

//6. Capital accumulation equation
  exp(li) = exp(lk)-(1-delta)*exp(lk(-1));


//7. Capital rental rate
 exp(Rk) = alpha*exp(ly)/exp(lk(-1));


//8. TFP shock
  z = rho*z(-1)+e;

end;

steady_state_model;
N = 1/3;
KoverN = (alpha/(1/beta-1+delta))^(1/(1-alpha));
Y = (KoverN)^alpha*N;
I = delta*KoverN*N;
K = KoverN*N;
C = Y-I;
w = (1-alpha)*K^alpha*N^(-alpha);
Rk1 = alpha*Y/K;
chi=(1-alpha)*(KoverN)^alpha*(1-N)/C; %for KPR0 and KPR3

ly=log(Y); lc=log(C); lk=log(K);
lh=log(N); li=log(I); lw=log(w);
Rk=log(Rk1); 
z=0;

end;
 
//steady;



%----------------------------------------------------------------
% 5. Estimation
%----------------------------------------------------------------


estimated_params;
alpha, beta_pdf, 0.35, 0.02;
beta, beta_pdf, 0.99, 0.002;
delta, beta_pdf, 0.025, 0.003;
rho, beta_pdf, 0.9, 0.05;
stderr e, inv_gamma_pdf, 0.01, inf;
end;

estimation(datafile=simuldataRBC,nobs=400,order=1,first_obs=500,mh_replic=2000,mh_nblocks=2,mh_drop=0.45,mh_jscale=0.8,mode_compute=6) ly;
%stoch_simul(periods=1000, order = 1);
