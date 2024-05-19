%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RBC model with stochastic trends
% Technology Shocks with Unit Root
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc


var k c N i y w Rk a;
varexo e;

parameters alpha beta delta chi muz sigmae;

% Calibration
alpha = 0.33;
beta = 0.99;
muz = 0.0033;
delta = 0.025;


[param,ss] = calibtrend;
chi = param(4);

sigmae = 0.01;

model;
//Consumption Euler equation
1/exp(c) = beta/exp(c(+1))*exp(a(+1)/(alpha-1))*(exp(Rk(+1))+1-delta);

//Labor supply
chi*exp(c)/(1-exp(N)) = exp(w);

//Labor demand
exp(w) = (1-alpha)*exp(alpha/(alpha-1)*a)*exp(k(-1)*alpha)*exp(N*(-alpha));

//Output
exp(y) = (exp(a)^(-alpha/(1-alpha)))*exp(k(-1))^(alpha)*exp(N)^(1-alpha);

//Law of motion for capital
exp(k) = exp(1/(alpha-1)*a)*(1-delta)*exp(k(-1))+ exp(i);

//Rental rate
exp(Rk) = alpha*exp(y)/exp(k(-1))*exp(a/(1-alpha));

//Resource constraint
exp(i) = exp(y) - exp(c);

//Shocks
a = muz + e;
end;

initval;
k = log(ss(1));
c = log(ss(2));
N = log(ss(3));
i = log(ss(4));
y = log(ss(5));
w = log(ss(6));
Rk = log(ss(7));
a = 0.0033;
end;

shocks;
var e = sigmae^2;
end;

steady;

H = 40;
stoch_simul(order = 1, irf = 40);

a_e = cumsum(a_e);
c_e = c_e + a_e*(1/(1-alpha));
y_e = y_e + a_e*(1/(1-alpha));
k_e = k_e + a_e*(1/(1-alpha));
i_e = i_e + a_e*(1/(1-alpha));
w_e = w_e + a_e*(1/(1-alpha));

k_e2(1,1) = 0;

for j = 2:H
k_e2(j,1) = k_e(j-1,1);
end


time=1:1:40;

figure
subplot(3,2,1);
plot(time,y_e, 'LineWidth',2)
title('Y');

subplot(3,2,2);
plot(time,c_e,'LineWidth',2)
title('C');

subplot(3,2,3)
plot(time,i_e,'LineWidth',2)
title('I');

subplot(3,2,4)
plot(time, N_e,'LineWidth',2)
title('N')

subplot(3,2,5)
plot(time, Rk_e,'LineWidth',2)
title('Rk')

subplot(3,2,6)
plot(time,w_e,'LineWidth',2)
title('w')

print -depsc figure_14_10.eps
