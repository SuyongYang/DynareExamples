// This code computes the basic DNK model in Chapter 19. 
// Plots Figures 19.1 and 19.2.
// @by Jianjun Miao

close all
var x pi i r_n v rr;
varexo e_r e_v;
parameters beta gamma nu theta phi_y phi_pi rho_r rho_v kappa ;

beta =0.99;
gamma = 1;
nu =1;
theta = 0.8;
phi_y = 0.25;
phi_pi = 1.5;
rho_r = 0.5;
rho_v = 0.5;

kappa =(1 - theta )*(1 - theta * beta )/ theta *(gamma + nu );

model (linear); //No need to let Dynare linearize
x=x (+1) - 1/gamma *(i-pi (+1) - r_n ); // IS curve
pi= kappa *x+ beta *pi (+1);   //NKPC curve
i= phi_pi *pi+ phi_y *x +v; //Taylor rule
rr=i-pi(+1); // Real interest rate
r_n = rho_r *r_n ( -1)+ e_r; // Natural rate shock
v = rho_v * v(-1)+e_v; // Monetary policy shock
end ;

shocks ;
var e_r =0.01^2;
var e_v = 0.01^2;
end ;

stoch_simul ( irf =11);

close all

time = 0:1:10

figure //Plot impulse responses to monetary policy shock
plot(time, x_e_v*100, time, pi_e_v*100, '--', time, i_e_v*100, '-.', time, rr_e_v*100,':','LineWidth',1.8)
legend('Output gap','Inflation rate','Nominal interest rate','Real interest rate')
xlabel('Time')
ylabel('Percent')
hold on
plot(time,zeros(11,1))
print -depsc2 fig19.1.eps

figure //Plot impulse responses to natural rate shock
plot(time, x_e_r*100, time, pi_e_r*100, '--', time, i_e_r*100, '-.', time, rr_e_r*100,':','LineWidth',1.8)
legend('Output gap','Inflation rate','Nominal interest rate','Real interest rate')
xlabel('Time')
ylabel('Percent')
hold on
plot(time,zeros(11,1))
print -depsc2 fig19.2.eps

