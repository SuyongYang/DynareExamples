% Basic RBC Model
% gauthier[at]vermandel.fr

close all;
%format long

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c h w r a g;
varexo e_a e_g;

parameters	beta alpha delta sigmaC sigmaL gy hc
			rho_a rho_g
			R Y C Z
			;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

% Loading calibration
@#include "calibration.mod"

%----------------------------------------------------------------
% 3. Model (the number refers to the equation in the paper)
%----------------------------------------------------------------
model(linear);
    %% Household
	% Euler
	sigmaC*(c(+1)-c)=r;
	% hours supply
	w = (1/sigmaL)*h+sigmaC*c;

    % Intermediary firms
	% Production function
	y = a + (1-alpha)*h;
	% Inputs Cost minimization
	w = y - h;
    % Resources constraint
	y = (1-gy)*c + gy*g ;
	
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
end;


stoch_simul(order = 1,irf=30) y c w h r;
epsilon_mat=zeros(options_.irf,size(var_list_,1)*size(M_.exo_names,1));
count = 0;
for i1=1:size(M_.exo_names,1)
	for i2=1:size(var_list_,1)
		count = count+1;
		epsilon_mat(:,count) = eval(['oo_.irfs.' var_list_(i2,:) '_' M_.exo_names(i1,:) '*100']);
	end
end
save(['irf_' M_.fname],'epsilon_mat');


