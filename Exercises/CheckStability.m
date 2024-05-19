%eigen value 

% Gap(i+1) = Gap(i) + (istar(i,1) - r(i,1) + (phi -1)*Infl(i))*dt;
% Infl(i+1) = Infl(i) + (rho*Infl(i) - kapa*Gap(i))*dt;

% x=[Gap Infl]'
% x(t+1) = [    1        (phi-1)dt  ] x(t)  +  [(istar-r)dt ]
%          [-kapa*dt      rho*dt+1  ]          [     0      ]

clear all

phi = 1.2;
dt = 0.01;
rho = 1/1.05;
ep = 1.05;
pai = 3.7;
theta =1/1.2;
kapa = (ep -1)*(1+pai)/theta;

% tem = rho - 4*kapa*(phi -1)
% tem01 = .5*(rho + tem^.5)*dt
% tem02 = .5*(rho - tem^.5)*dt

A = [1 (phi-1)*dt; - kapa*dt rho*dt+1]

eig(A)