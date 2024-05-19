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


%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var ly lc lk li lh lw Rk z;
varexo e;

parameters beta delta chi alpha rho;

%----------------------------------------------------------------
% 2. Declaring observable variables
%---------------------------------------------------------------

varobs ly;


%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
  
  //Consumption Euler equation
  (1/exp(lc)) = beta*(1/exp(lc(+1)))*(1+alpha*(exp(lk)^(alpha-1))*exp(z(+1))*exp(lh(+1))^(1-alpha)-delta);

// Labor supply
 chi*exp(lc)/(1-exp(lh)) = exp(lw);

// Labor demand
 exp(lw) = exp(z)*(1-alpha)*exp(lk(-1))^alpha*exp(lh)^(-alpha);

//Resource constraint
  exp(lc)+exp(li) = exp(ly);

//Production function
  exp(ly) = exp(z)*(exp(lk(-1))^alpha)*(exp(lh))^(1-alpha);

//Capital accumulation equation
  exp(li) = exp(lk)-(1-delta)*exp(lk(-1));


//Capital rental rate
 exp(Rk) = alpha*exp(ly)/exp(lk(-1));


//TFP shock
  z = rho*z(-1)+e;

end;

%-------------------------------------------------------------
% 4. Specifying Steady State
%--------------------------------------------------------------


initval;
  lk = 2.2;
  lc = -0.26;
  lh = -1.098;  
  li = -1.44;
  ly = 0.005;
  lw= 0.7;
  Rk= -3.3;
  z = 0; 
end;




%----------------------------------------------------------------
% 5. Estimation
%----------------------------------------------------------------


estimated_params;
alpha, beta_pdf, 0.35, 0.02;
beta, beta_pdf, 0.99, 0.002;
delta, beta_pdf, 0.025, 0.003;
chi, gamma_pdf, 1.65, 0.02;
rho, beta_pdf, 0.9, 0.05;
stderr e, inv_gamma_pdf, 0.01, inf;
end;

estimation(datafile=simuldataRBC,nobs=200,order=1,first_obs=500,mh_replic=2000,mh_nblocks=2,mh_drop=0.45,mh_jscale=0.8,mode_compute=6) ly;
