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

clear
addpath('../Tools')

load simulated_data_dynare_RBC_without_SS_2000_periods

%use 2000 periods
Y0=[c(1:2000,1) ghat(1:2000,1)];

%chose starting values
rhoz    = 0.5;  
rhog    = 0.5;  
eps_z   = 0.1;%0.0068;
eps_g   = 0.1;
par_names={'rhoz','rhog','eps_z','eps_g'};
[rhoz, rhog, eps_z, eps_g]=par_transform_unbounded(rhoz,rhog,eps_z,eps_g);
x_start=[rhoz rhog eps_z eps_g];


%% start Optimizer for ML
% initialize optimizer
H0 = 1e-4*eye(length(x_start));
crit = 1e-7;
nit = 1000;

[fhat,xhat] = csminwel(@likelihood_Kalman,x_start,H0,[],crit,nit,Y0);
%retransform variables and display them
[rhoz, rhog, eps_z, eps_g]=par_retransform_bounded(xhat(1),xhat(2),xhat(3),xhat(4));
xhat_retrans=[rhoz, rhog, eps_z, eps_g];
for ii=1:length(xhat)
    fprintf('%10s: %4.3f\n',par_names{ii},xhat_retrans(ii))
end