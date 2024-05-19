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

clear all
addpath('../../Matlab_tools')
rand('state',1)
randn('state',1)

load simulated_data_dynare_RBC_without_SS_2000_periods
Y0=[c ghat];

rhoz    = 0.95;  
rhog    = 0.95;  
eps_z   = 0.0068;
eps_g   = 0.01;
par_names={'rhoz','rhog','eps_z','eps_g'};
[rhoz, rhog, eps_z, eps_g]=par_transform_unbounded(rhoz,rhog,eps_z,eps_g);
x_start=[rhoz rhog eps_z eps_g];


%% Find Posterior Mode
% initialize optimizer
H0 = 1e-4*eye(length(x_start));
crit = 1e-7;
nit = 1000;
verbose = 2;

[fhat,xhat,ghat,inv_Hessian] = csminwel(@posterior_Kalman,x_start,H0,[],crit,nit,Y0);

% inv_Hessian = Hhat\eye(4);
Sigma_chol_lower = cholcov(inv_Hessian)';

%% Bayesian Estimation
MH_draws=10000;
burnin = 1000;
scale_mh = 1;
accept=0;
draws=zeros(MH_draws,4);
draws(1,:)=xhat;
old_posterior = -1*posterior_Kalman(xhat,Y0);
proposal_draws = scale_mh*Sigma_chol_lower*randn(4,MH_draws);% Le*(Le)'=L*L'=Omega
for ii=2:MH_draws
    xhatstar = draws(ii-1,:)+proposal_draws(:,ii)';
    new_posterior = -1*posterior_Kalman(xhatstar,Y0);
    accprob=exp(new_posterior-old_posterior);
    if rand(1,1)<=accprob
        draws(ii,:)=xhatstar;
        old_posterior = new_posterior;
        accept=accept+1;
    else
        draws(ii,:)=draws(ii-1,:);
    end
    if mod(ii,50)==0
        ratio=accept/ii;
        disp(['Acceptance Rate:' num2str(ratio)])
    end
end
[rhoz, rhog, eps_z, eps_g]=par_retransform_bounded(draws(burnin+1:end,1),draws(burnin+1:end,2),draws(burnin+1:end,3),draws(burnin+1:end,4));

figure
subplot(2,2,1)
hist(rhoz,25)
title('\rho_z')
subplot(2,2,2)
hist(rhog,25)
title('\rho_g')
subplot(2,2,3)
hist(eps_z,25)
title('\epsilon_z')
subplot(2,2,4)
hist(eps_g,25)
title('\epsilon_g')

[rhoz, rhog, eps_z, eps_g]=par_retransform_bounded((draws(burnin+1:end,1)),mean(draws(burnin+1:end,2)),mean(draws(burnin+1:end,3)),mean(draws(burnin+1:end,4)));