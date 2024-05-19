function [F,H,alarm_flag]=get_solution_RBC(rhoz,rhog)
% function [F,H,alarm_flag]=get_solution_RBC(rhoz,rhog)
% Returns the state space solution matrices F and H for the solved RBC
% model
%
% Inputs:
%   - rhoz      [scalar]                autocorrelation coefficient TFP
%   - rhog      [scalar]                autocorrelation coefficient G
%
% Outputs:
%   - F         [n_state by n_state]    state transition matrix
%   - H         [n_control by n_state]  matrix mapping controls into states
%   - alarm_flag [scalar]               1= Blanchard-Kahn conditions not
%                                           satisfied

% Copyright (C) 2014-2018 Benjamin Born and Johannes Pfeifer
% 
%  This is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  It is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  For a copy of the GNU General Public License,
%  see <http://www.gnu.org/licenses/>.

%% RBC_Model, Calibration for US
%% define parameters and do model calibration
sigma=1;
alpha= 0.33;
i_y=0.25;
k_y=10.4;
x=0.0055;
n=0.0027;
gshare= 0.2038;
gammax=(1+n)*(1+x);
delta=i_y/k_y-x-n-n*x;
beta=(1+x)*(1+n)/(alpha/k_y+(1-delta));

%% compute steady states
l_ss=0.33;
k_ss = ((1/beta*(1+n)*(1+x)-(1-delta))/alpha)^(1/(alpha-1))*l_ss; 
i_ss = (x+n+delta+n*x)*k_ss;
y_ss=k_ss^alpha*l_ss^(1-alpha);
g_ss=gshare*y_ss;
c_ss = (1-gshare)*k_ss^(alpha)*l_ss^(1-alpha)-i_ss;
psi=(1-alpha)*(k_ss/l_ss)^alpha*(1-l_ss)/c_ss^sigma;

%% define coefficients and auxiliary variables
gamma_l=(l_ss/(1-l_ss)+alpha)^-1;

alpha1=1/gammax*(alpha*y_ss/k_ss*(1+(1-alpha)*gamma_l)+(1-delta));
alpha2=-(c_ss/(gammax*k_ss)+y_ss/(gammax*k_ss)*(1-alpha)*gamma_l*sigma);
alpha3=y_ss/(gammax*k_ss)*(1+(1-alpha)*gamma_l);
alpha4=-g_ss/(gammax*k_ss);

alpha5=-beta/(gammax*sigma)*alpha*(k_ss/l_ss)^(alpha-1)*(alpha-1)*(1-alpha*gamma_l);
alpha6=1-beta/gammax*alpha*(k_ss/l_ss)^(alpha-1)*(alpha-1)*gamma_l;
alpha7=-beta/(gammax*sigma)*alpha*(k_ss/l_ss)^(alpha-1)*(1-(alpha-1)*gamma_l);

%% Solve Model using eigenvalue decomposition
%% matrices of linear RE model: A*E_t(x_{t+1})=B*(x_{t})
% here; x=[k_t z_t g_t c_t]
A= [1       0           0           0; 
    alpha5  alpha7      0           alpha6;
    0       0           1           0 
    0       1           0           0
   ];
    
B= [alpha1  alpha3  alpha4  alpha2 ;
    0       0           0       1
    0       0           rhog    0;
    0       rhoz        0       0];


%% multiplying with inv(A) 
W=A\B;

%% partition matricies (the first three elements in x are states)
%% we need w11 and w12 for M
n_x=3;
n_x1=1; %endogenous states
n_x2=2; %exogenous states
n_u=1;

w11=W(1:n_x,1:n_x);
w12=W(1:n_x,n_x+1:end);

%% jordan decomposition: W=D*\Lambda*D^{-1}
[evec,eval]=eig(W);
[mu,ind]=sort(abs(diag(eval)));% sort eigenvalues in increasing order
dinv=inv(evec(:,ind));         % take inverse of matrix D

%% partition D^{-1} 
q21=dinv(n_x+1:end,1:n_x);               
q22=dinv(n_x+1:end,n_x+1:end);                 

%% Model solution 
F = w11-w12*q21/q22;            % transition under optimal policy 
H = -q22\q21;      % policy function mapping state into control

if any(any(isnan([F H']))) || any(any(isinf([F H'])))
    alarm_flag=1;
else
    alarm_flag=0;
end