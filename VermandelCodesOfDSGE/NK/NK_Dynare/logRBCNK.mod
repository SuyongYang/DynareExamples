% Basic RBC Model
% With New Keynesian Features
% + Sluggish prices (forward)
% Vermandel 2013
% gauthier[at]vermandel.fr

close all;
%format long

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y
	c
	r
	pi
	w z
	h k i
	mc
	a g
	;

varexo	e_a e_g e_r;

parameters	beta alpha delta sigmaC sigmaL gy
			rho_r phi_r phi_y theta_p
			rho_a rho_g R Y Z C I
			
			;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
@#include "calibration.mod"


%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model(linear);
    %% Household
	% Euler
	sigmaC*(c(+1)-c)=r-pi(+1);
	% hours supply
	w=(1/sigmaL)*h+sigmaC*c;

	%% Capital supply
	% No arbitrage Bonds-Capital
	R*(r-pi(+1)) = Z*z(+1);
	% Capital law of motion
	delta*i = k-(1-delta)*k(-1);

    % Intermediary firms
	% Production function
	y = a + alpha*k(-1) + (1-alpha)*h;
	% Real marginal cost
	mc = alpha*z + (1-alpha)*w - a;
	% Cost minimization
	w + h = z + k(-1);
	% Price dynamics NKPC
	pi = beta*pi(+1) + (1-theta_p)*(1-theta_p*beta)/theta_p*mc;
	
    % Resources constraint
	Y*y = C*c + I*i + gy*Y*g ;
	
	% monetary policy
	r = rho_r*r(-1) + (1-rho_r)*(phi_r*pi +  phi_y*(y-y(-1))) + e_r;
	
    % Exogenous shocks
	a = rho_a*a(-1) + e_a;
	g = rho_g*g(-1) + e_g;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
check;
steady;

shocks;
var e_a;  stderr .01;
var e_g;  stderr .01;
var e_r;  stderr .01;
end;

stoch_simul(order=1,irf=50);