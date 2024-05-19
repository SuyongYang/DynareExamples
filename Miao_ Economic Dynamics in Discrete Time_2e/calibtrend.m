%
% Calibration of the basic RBC model with KPR0 and KPR3 Preferences
%

function [param,SS] = calibtrend

alpha=0.33;
beta=0.99;
delta =0.025;
N=0.33;
muz = 0.0033;

Rk = 1/beta/exp(muz/(alpha-1))-1+delta;

KoverN = (alpha/Rk/exp(muz/(1-alpha)))^(1/(1-alpha));
N=1/3;
Y = (KoverN)^alpha*N*exp(alpha/(alpha-1)*muz);
I = KoverN*N*(1-exp(muz/(alpha-1))*(1-delta));
K = KoverN*N;
C = Y-I;

w = (1-alpha)*K^alpha*N^(-alpha)*exp(alpha/(alpha-1)*muz);
chi=w*(1-N)/C; %for KPR0 and KPR3

param = [alpha, beta, delta, chi];
SS = [K,C,N,I,Y,w, Rk];
