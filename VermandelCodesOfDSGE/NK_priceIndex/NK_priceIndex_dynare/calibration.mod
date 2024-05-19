
alpha   	= 0.36;				% share of capital in ouput
beta    	= 0.99;				% discount factor
delta  		= 0.025;			% depreciation of capital
sigmaC		= 1;				% risk aversion consumption
sigmaL		= 2;				% labor disutility
gy 			= 0.2; 			 	% Public spending in GDP
theta_p		= .75;				% new keynesian Philips Curve, forward term
gamma_p		= .2;				% new keynesian Philips Curve, backward term
epsilon_p	= 10;				% subsituability/mark-up on prices
rho_r		= .7;				% Monetary Policy Smoothing Parameter
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
