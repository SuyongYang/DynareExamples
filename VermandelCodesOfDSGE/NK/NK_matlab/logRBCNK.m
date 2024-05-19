% Basic RBC Model
% With New Keynesian Features
% + Sluggish prices (forward)
% gauthier[at]vermandel.fr

clear all
close all 

T = 50; % Number of periods for the IRF

%%%%%%%%%%%%%%%%%%%% CALIBRATION %%%%%%%%%%%%%%%%%%%%
% RBC Parameters
alpha   	= 0.36;				% share of capital in ouput
beta    	= 0.99;				% discount factor
delta  		= 0.025;			% depreciation of capital
sigmaC		= 1;				% risk aversion consumption
sigmaL		= 2;				% labor disutility
gy 			= 0.2; 			 	% Public spending to GDP

% New Keynesian Parameters
theta_p		= .75;				% new keynesian Philips Curve, forward term
epsilon_p	= 10;				% subsituability/mark-up on prices
rho_r		= 0.7;				% Monetary Policy Smoothing Parameter
phi_y		= .125;				% Monetary Policy GDP Growth Target
phi_r		= 1.5;				% Monetary Policy Inflation Growth Target

% shock process
rho_a   	= 0.95; 			% productivity 
rho_g   	= 0.95; 			% public spending

% steady states
R			= 1/beta;
Z			= 1/beta-(1-delta);
H			= 1/3;				
MC			= (epsilon_p-1)/epsilon_p;
K			= H*(Z/(alpha*MC))^(1/(alpha-1));
Y			= K^alpha*H^(1-alpha);
C			= (1-gy)*Y-delta*K;
I			= delta*K;
W			= (1-alpha)*MC*Y/H;
DELTAC		= C^-sigmaC;
chi			= DELTAC*W/(H^(1/sigmaL));

%%%%%%%%%%%%%%%%%%%% Variable ID assignments %%%%%%%%%%%%%%%%%%%%
% Standard variables X(t)
var_y	= 1;		% Production
var_c	= 2;		% Consumption
var_i	= 3;		% Investment
var_h	= 4;		% Hours
var_k	= 5;		% Physical Capital
var_w	= 6;		% Real Wage
var_z	= 7;		% Real Capital Cost
var_r	= 8;		% Real Rate
var_pi	= 9;		% Inflation Rate
var_mc	= 10;		% Real Marginal Cost
var_s_a	= 11;		% Productivity Shock Process
var_s_g	= 12;		% Public Spending Shock Process
% Lagged variables X(t-1)
var_s_a_lag	= 13;	% Productivity Shock Process
var_s_g_lag	= 14;	% Public Spending Shock Process
var_r_lag	= 15;	% Lagged Rate
var_k_lag	= 16;	% Lagged Physical Capital
var_y_lag	= 17;	% Lagged Production

% Shock ID assignment
e_a	= 1;	% Productivity Shock
e_g	= 2;	% Spending Shock
e_r	= 3;	% Monetary Policy Shock

%%%%%%%%%%%%%%%%%%%% MATRIX DEFINITIONS %%%%%%%%%%%%%%%%%%%%
% number of endogenous variables
nx = 17; 
% Number of shocks
nz = 3; 
nu = nz;
% number of predetermined variables
nk = 5;
% (do not edit) : declaring matrix
AA = zeros(nx,nx);
BB = zeros(nx,nx);
CC = zeros(nx,nu);
VC = zeros(nu,nu);

%%%%%%%%%%%%%%%%%%%% VARIANCE-COVARIANCE MATRIX %%%%%%%%%%%%%%%%%%%%
VC(1,1) = 0.01^2; % 1% productivity shock
VC(2,2) = 0.01^2; % 1% spending shock
VC(3,3) = 0.01^2; % 1% rate shock


%%%%%%%%%%%%%%%%%%%% MODEL EQUATIONS %%%%%%%%%%%%%%%%%%%%
EqNb 				= 0; % starting the counter;

% Euler
EqNb 				= EqNb+1; % equation number
AA(EqNb,var_c)		= sigmaC;
BB(EqNb,var_c)		= -sigmaC;
BB(EqNb,var_r)		= -1;
AA(EqNb,var_pi)		= 1;

% Hours supply
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_w)		= -1;
BB(EqNb,var_h)		= 1/sigmaL;
BB(EqNb,var_c)		= sigmaC;

% No Arbitrage Capital Bonds
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_r)		= R;
AA(EqNb,var_pi)		= -R;
AA(EqNb,var_z)		= -Z;

% Capital Law of Motion
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_i)		= -delta;
BB(EqNb,var_k)		= 1;
BB(EqNb,var_k_lag)	= -(1-delta);

% Production function
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_y)		= -1;
BB(EqNb,var_s_a)	= 1;
BB(EqNb,var_k_lag)	= alpha;
BB(EqNb,var_h)		= 1-alpha;

% Marginal Cost of Production
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_mc)		= -1;
BB(EqNb,var_z)		= alpha;
BB(EqNb,var_w)		= 1-alpha;
BB(EqNb,var_s_a)	= -1;

% Cost minimization
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_w)		= -1;
BB(EqNb,var_h)		= -1;
BB(EqNb,var_k_lag)	= 1;
BB(EqNb,var_z)		= 1;

% Price dynamics NKPC
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_pi)		= -1;
AA(EqNb,var_pi)		= beta;
BB(EqNb,var_mc)		= (1-theta_p)*(1-theta_p*beta)/theta_p;
	
% Resources constraint
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_y)		= -Y;
BB(EqNb,var_c)		= C;
BB(EqNb,var_i)		= I;
BB(EqNb,var_s_g)	= gy*Y;
	
% Monetary Policy
EqNb 					= EqNb+1; % equation number
BB(EqNb,var_r)		= -1;
BB(EqNb,var_r_lag)	= rho_r;
BB(EqNb,var_pi)		= (1-rho_r)*phi_r;
BB(EqNb,var_y)		= (1-rho_r)*phi_y;
BB(EqNb,var_y_lag)	= -(1-rho_r)*phi_y;
CC(EqNb,e_r)		= 1;

% Shock process : Productivity
EqNb 					= EqNb+1; % equation number
BB(EqNb,var_s_a)		= -1;
BB(EqNb,var_s_a_lag)	= rho_a;
CC(EqNb,e_a)			= 1;

% Shock process : Spending
EqNb 					= EqNb+1; % equation number
BB(EqNb,var_s_g)		= -1;
BB(EqNb,var_s_g_lag)	= rho_g;
CC(EqNb,e_g)			= 1;


% DEFINING LAGGED VARIABLES
%Define s_a(t-1)
EqNb 					= EqNb+1; % equation number
AA(EqNb,var_s_a_lag)	= -1;
BB(EqNb,var_s_a)		= 1;

%Define s_g(t-1)
EqNb 					= EqNb+1; % equation number
AA(EqNb,var_s_g_lag) 	= -1;
BB(EqNb,var_s_g) 		= 1;

%Define R(t-1)
EqNb 					= EqNb+1; % equation number
AA(EqNb,var_r_lag)	 	= -1;
BB(EqNb,var_r)	 		= 1;

%Define K(t-1)
EqNb 					= EqNb+1; % equation number
AA(EqNb,var_k_lag)	 	= -1;
BB(EqNb,var_k)	 		= 1;

%Define Y(t-1)
EqNb 					= EqNb+1; % equation number
AA(EqNb,var_y_lag)	 	= -1;
BB(EqNb,var_y)	 		= 1;

% Solving the model
[m1,m2,m3,m4]=solvek(-AA,BB,CC,VC,nk);

% Plotting results
for i1 = 1:nz
	mprea=zeros(nk,1);
	ma=zeros(nx-nk,1);
	exog=zeros(nz,1);
	exog(i1) = 1*sqrt(VC(i1,i1));
	for i=2:(T+1)
		ma(:,i)=m1*mprea(:,i-1)+m2*exog;
		mprea(:,i)=m3*mprea(:,i-1)+m4*exog;
		exog(i1) = 0;
	end

	figure
	subplot(3,3,1)
	plot(1:T,ma(var_y,2:(T+1)))
	xlim([1 T])
	title('y')
	
	subplot(3,3,2)
	plot(1:T,ma(var_c,2:(T+1)))
	xlim([1 T])
	title('c')
	
	subplot(3,3,3)
	plot(1:T,ma(var_i,2:(T+1)))
	xlim([1 T])
	title('i')
	
	subplot(3,3,4)
	plot(1:T,ma(var_r,2:(T+1)))
	xlim([1 T])
	title('r')
	
	subplot(3,3,5)
	plot(1:T,ma(var_h,2:(T+1)))
	xlim([1 T])
	title('h')
	
	subplot(3,3,6)
	plot(1:T,ma(var_w,2:(T+1)))
	xlim([1 30])
	title('w')
	
	subplot(3,3,7)
	plot(1:T,ma(var_pi,2:(T+1)))
	xlim([1 T])
	title('\pi')

end






