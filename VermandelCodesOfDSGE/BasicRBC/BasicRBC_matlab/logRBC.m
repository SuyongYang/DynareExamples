clear all
close all 
%%%%%%%%%%%%%%%%%%%% CALIBRATION %%%%%%%%%%%%%%%%%%%%
% RBC Parameters
beta 	= 0.99; 	%Discount Factor
alpha 	= 0.36;		%Capital share
gy 		= 0.2;    	%Public spending in GDP
sigmaC 	= 2.5;		%Risk aversion
sigmaL 	= .5; 		%Elasticity of labor

% autoregressive roots parameters
rho_a	= 0.95;
rho_g	= 0.95;


%%%%%%%%%%%%%%%%%%%% Variable ID assignments %%%%%%%%%%%%%%%%%%%%
% Standard variables X(t)
var_y	= 1;		% Production
var_c	= 2;		% Consumption
var_h	= 3;		% Hours
var_w	= 4;		% Real Wage
var_r	= 5;		% Real Rate
var_s_a	= 6;		% Productivity Shock Process
var_s_g	= 7;		% Public Spending Shock Process
% Lagged variables X(t-1)
var_s_a_lag	= 8;	% Productivity Shock Process
var_s_g_lag	= 9;	% Public Spending Shock Process

% Shock ID assignment
e_a	= 1;	% Productivity Shock
e_g	= 2;	% Spending Shock

%%%%%%%%%%%%%%%%%%%% MATRIX DEFINITIONS %%%%%%%%%%%%%%%%%%%%
% number of endogenous variables
nx = 9; 
% Number of shocks
nz = 2; 
nu = nz;
% number of predetermined variables
nk = 2;
% (do not edit) : declaring matrix
AA = zeros(nx,nx);
BB = zeros(nx,nx);
CC = zeros(nx,nu);
VC = zeros(nu,nu);

%%%%%%%%%%%%%%%%%%%% VARIANCE-COVARIANCE MATRIX %%%%%%%%%%%%%%%%%%%%
VC(1,1) = 0.01^2; % 1% productivity shock
VC(2,2) = 0.01^2; % 1% spending shock


%%%%%%%%%%%%%%%%%%%% MODEL EQUATIONS %%%%%%%%%%%%%%%%%%%%
EqNb 				= 0; % starting the counter;

% Euler
EqNb 				= EqNb+1; % equation number
AA(EqNb,var_c)		= sigmaC;
BB(EqNb,var_c)		= -sigmaC;
BB(EqNb,var_r)		= -1;

% Hours supply
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_w)		= -1;
BB(EqNb,var_h)		= 1/sigmaL;
BB(EqNb,var_c)		= sigmaC;

% Production function
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_y)		= -1;
BB(EqNb,var_s_a)	= 1;
BB(EqNb,var_h)		= 1-alpha;

% Cost minimization
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_w)		= -1;
BB(EqNb,var_y)		= 1;
BB(EqNb,var_h)		= -1;

% Resources constraint
EqNb 				= EqNb+1; % equation number
BB(EqNb,var_y)		= -1;
BB(EqNb,var_c)		= (1-gy);
BB(EqNb,var_s_g)	= gy;

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


% Solving the model
[m1,m2,m3,m4]=solvek(-AA,BB,CC,VC,nk);

% Plotting results
for i1 = 1:nz
	mprea=zeros(nk,1);
	ma=zeros(nx-nk,1);
	exog=zeros(nz,1);
	exog(i1) = 1*sqrt(VC(i1,i1));
	for i=2:50
		ma(:,i)=m1*mprea(:,i-1)+m2*exog;
		mprea(:,i)=m3*mprea(:,i-1)+m4*exog;
		exog(i1) = 0;
	end

	figure
	subplot(3,2,1)
	plot(1:30,ma(var_y,2:31))
	title('y')
	subplot(3,2,2)
	plot(1:30,ma(var_c,2:31))
	title('c')
	subplot(3,2,3)
	plot(1:30,ma(var_r,2:31))
	title('r')
	subplot(3,2,4)
	plot(1:30,ma(var_h,2:31))
	title('h')
	subplot(3,2,5)
	plot(1:30,ma(var_w,2:31))
	title('w')

end






