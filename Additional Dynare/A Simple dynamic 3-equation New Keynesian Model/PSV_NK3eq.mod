% The analytics of the New Keynesian 3-equation Model (2015)
% code by Gauthier Vermandel



%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y pi r s_D s_S s_r;

varexo e_D e_S e_R;

parameters  beta sigma phi chi rho phi_pi phi_y theta rho_D rho_S rho_R;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

beta    	= 0.99;				% discount factor
sigma		= 1;				% risk aversion consumption
phi			= 1;				% labor disutility
theta		= 3/4;				% new keynesian Philips Curve, forward term
epsilon		= 6;				% subsituability/mark-up on prices
rho			= 0;				% MPR Smoothing
phi_pi		= 1.5;				% MPR Inflation
phi_y		= 0.5/4;			% MPR GDP

% shock processes
rho_D   = 0.9;
rho_S   = 0.9;
rho_R 	= 0.4;

% steady states
R		= 1/beta;
H		= 1/3;
MC		= (epsilon-1)/epsilon;
W		= MC;
Y		= H;
chi		= W*Y^-sigma*H^-phi;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model(linear); 
	% IS curve
	y = y(+1) - 1/sigma*(r-pi(+1)) + s_D;
	% AS curve
	pi = beta*pi(+1) + ((1-theta)*(1-beta*theta)/theta)*(sigma+phi)*y + s_S;
	% Monetary Policy Rule
	r = rho*r(-1) +  (1-rho)*( phi_pi*pi + phi_y*y ) + s_r  ;
	
	% shocks
	s_D = rho_D*s_D(-1)+e_D;
	s_S = rho_S*s_S(-1)+e_S;
	s_r = rho_R*s_r(-1)+e_R;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
shocks;
var e_D;  stderr 1;
var e_S;  stderr 1;
var e_R;  stderr 1;
end;

%check;

% GENERATING IRF
% benchmark
set_param_value('phi_pi',1.5);
set_param_value('phi_y',0.5/4);
stoch_simul(irf=20,nograph) y r pi;
irf1 = oo_.irfs;
% inflation target
set_param_value('phi_pi',1.7);
set_param_value('phi_y',0.5/4);
stoch_simul(irf=20,nograph) y r pi;
irf2 = oo_.irfs;
% output target
set_param_value('phi_pi',1.5);
set_param_value('phi_y',0.8/4);
stoch_simul(irf=20,nograph) y r pi;
irf3 = oo_.irfs;

% near-flexible price
%set_param_value('theta',0.01);
%set_param_value('phi_y',1);
%stoch_simul(order=1,irf=30) y r pi;

figure;
set(gcf,'numbertitle','off','name','Demand shock','Position', [450, 450, 1000,200]) 
subplot(1,3,1)
plot(1:options_.irf,irf1.y_e_D,'gx-')
hold on;
plot(1:options_.irf,irf2.y_e_D,'ro-')
plot(1:options_.irf,irf3.y_e_D,'b--')
hold off;
title('y_t - output gap')
grid on;
subplot(1,3,2)
plot(1:options_.irf,irf1.pi_e_D,'gx-')
hold on;
plot(1:options_.irf,irf2.pi_e_D,'ro-')
plot(1:options_.irf,irf3.pi_e_D,'b--')
hold off;
title('\pi_t - inflation')
grid on;
subplot(1,3,3)
plot(1:options_.irf,irf1.r_e_D,'gx-')
hold on;
plot(1:options_.irf,irf2.r_e_D,'ro-')
plot(1:options_.irf,irf3.r_e_D,'b--')
hold off;
title('r_t - interest rate')
grid on;

figure;
set(gcf,'numbertitle','off','name','Supply shock','Position', [450, 450, 1000,200]) 
subplot(1,3,1)
plot(1:options_.irf,irf1.y_e_S,'gx-')
hold on;
plot(1:options_.irf,irf2.y_e_S,'ro-')
plot(1:options_.irf,irf3.y_e_S,'b--')
hold off;
title('y_t - output gap')
grid on;
subplot(1,3,2)
plot(1:options_.irf,irf1.pi_e_S,'gx-')
hold on;
plot(1:options_.irf,irf2.pi_e_S,'ro-')
plot(1:options_.irf,irf3.pi_e_S,'b--')
hold off;
title('\pi_t - inflation')
grid on;
subplot(1,3,3)
plot(1:options_.irf,irf1.r_e_S,'gx-')
hold on;
plot(1:options_.irf,irf2.r_e_S,'ro-')
plot(1:options_.irf,irf3.r_e_S,'b--')
hold off;
title('r_t - interest rate')
grid on;


figure;
set(gcf,'numbertitle','off','name','Monetary Policy shock','Position', [450, 450, 1000,200]) 
subplot(1,3,1)
plot(1:options_.irf,irf1.y_e_R,'gx-')
hold on;
plot(1:options_.irf,irf2.y_e_R,'ro-')
plot(1:options_.irf,irf3.y_e_R,'b--')
hold off;
title('y_t - output gap')
grid on;
subplot(1,3,2)
plot(1:options_.irf,irf1.pi_e_R,'gx-')
hold on;
plot(1:options_.irf,irf2.pi_e_R,'ro-')
plot(1:options_.irf,irf3.pi_e_R,'b--')
hold off;
title('\pi_t - inflation')
grid on;
subplot(1,3,3)
plot(1:options_.irf,irf1.r_e_R,'gx-')
hold on;
plot(1:options_.irf,irf2.r_e_R,'ro-')
plot(1:options_.irf,irf3.r_e_R,'b--')
hold off;
title('r_t - interest rate')
grid on;




