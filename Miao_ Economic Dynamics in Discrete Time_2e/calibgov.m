%
% Calibration of the basic RBC model with Government Spending 
%

function [param,SS] = calibgov

alpha=0.33;
beta=0.99;
delta =0.025;
N=0.33;

Rk = 1/beta-1+delta;
KoverN = (alpha/Rk)^(1/(1-alpha));
N=1/3;
Y = (KoverN)^alpha*N;
G=0.2*Y;
I = delta*KoverN*N;
K = KoverN*N;
C = Y-I-G;
chi=(1-alpha)*(KoverN)^alpha*(1-N)/C; %for KPR0 and KPR3

w = (1-alpha)*K^alpha*N^(-alpha);

Rs = Rk+1-delta
Rf = 1/beta; 

param = [alpha, beta, delta, chi];
SS = [K,C,N,I,Y,w, Rk,Rs,Rf,G];
