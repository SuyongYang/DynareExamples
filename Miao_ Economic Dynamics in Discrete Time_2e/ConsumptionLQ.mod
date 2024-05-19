% Basic RBC Model 
% Linearization in level
%
% 
%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;


%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var c a y z;
%var c a y;
varexo e ep;
%varexo e;

parameters beta r rho sigma sigmap;

%----------------------------------------------------------------
% 2. Quadratic Utility
%----------------------------------------------------------------
% 
% U(c,n) = -(c-b)^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load parameterfile
%set_param_value('rho',rho)

r      =0.04;
beta   = 1/(1+r);
sigma   = 0.01;
rho     = 0.8 ;
sigmap = 0.01;
%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
//(1) Euler equation
    c(+1)=c;
//(2) Budget constraint 
    c+a=(1+r)*a(-1)+y;
//(3) Income shock
   // y = rho*y(-1)+e;
 y = z + e;
 z = z(-1) + ep;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  
  c = 0;
  a = 0;
  y = 0; 
 z =0;
end;
steady;

check;

shocks;
var e = sigma^2;
var ep =sigmap^2;
end;


//stoch_simul(order=1, pruning,irf=40);
stoch_simul(order=1, irf=40);
