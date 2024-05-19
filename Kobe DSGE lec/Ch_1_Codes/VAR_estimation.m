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

clear all;

%% Data preparation

real_gdp = xlsread('BP.xlsx',1,'B8:B292');
gdp_deflator = xlsread('BP.xlsx',1,'D8:D292');
nom_gov_cons_inv = xlsread('BP.xlsx',1,'F8:F292');
nom_priv_cons_nd = xlsread('BP.xlsx',1,'J8:J292');
nom_priv_cons_serv = xlsread('BP.xlsx',1,'L8:L292');
nipa_pop=[NaN(4,1); xlsread('BP.xlsx','N8:N288')];

% construct real per capita values
Y = real_gdp./nipa_pop;
G = nom_gov_cons_inv./nipa_pop./gdp_deflator;
C = (nom_priv_cons_nd+nom_priv_cons_serv)./nipa_pop./gdp_deflator;

timeline = (1947:0.25:2018)';

x =[log(G) log(Y) log(C)];
var_names={'G','Y','C'};
G_pos=1;
Y_pos=2;
C_pos=3;

%restrict to sample before financial crisis
x=x(timeline>1953.75 & timeline<2008,:);
% x=x(timeline>=1948 & timeline<2018.25,:);


addpath('../../Matlab_tools/')

number_of_lags = 4;
constant_dummy =1;
trend_dummy = 1;
Z=[];

for lagnumber=1:number_of_lags
    tempmat = lagmatrix(x,lagnumber);
    Z = [Z tempmat(number_of_lags+1:end,:)];
end
[T, nvars] =size(x);

%add trend and constant
if trend_dummy==1
    Z = [(1:T-number_of_lags)' Z];
end

if constant_dummy==1
    Z = [ones(T-number_of_lags,1) Z];
end

%get OLS starting values
Z = Z';
Y=x(number_of_lags+1:end,:)';
T=size(Z,2);
bhat = kron((Z*Z')\Z,eye(nvars))*Y(:);
betta_OLS=reshape(bhat,nvars,nvars*number_of_lags+constant_dummy+trend_dummy);

resids=(eye(T)-Z'/(Z*Z')*Z)*Y';
L_OLS=chol(cov(resids),'lower');

AR_pars=(constant_dummy+trend_dummy)*nvars+nvars^2*number_of_lags;
total_pars=AR_pars+nvars*(nvars+1)/2;

x_start=zeros(total_pars,1);
x_start(1:AR_pars)=betta_OLS(:);
x_start(AR_pars+1:total_pars)=cholmat2vec(L_OLS);

H0=1e-4*eye(total_pars);
crit=1e-7;
nit=1000;

%run ML estimation
[fhat,xhat]=csminwel(@VAR_ML,x_start,H0,[],crit,nit,Y,Z);

%reconstruct parameters
betta=reshape(xhat(1:AR_pars,1),nvars,...
    nvars*number_of_lags+constant_dummy+trend_dummy);
L=vec2cholmat(xhat(AR_pars+1:end,1));
Sigma_U=L*L';

%construct IRFs
[F,Q]=build_companionform(betta(:,constant_dummy+trend_dummy+1:end),number_of_lags,...
    L);

n_augmented_vars=size(F,1);
IRF_length=20;
IRFs=zeros(n_augmented_vars,IRF_length);

shock_vector=zeros(n_augmented_vars,1);
shock_pos=G_pos;
shock_vector(shock_pos)=1/Q(shock_pos,shock_pos); %set to 1 percent of GDP

for ii=1:IRF_length
    IRFs(:,ii)=F^(ii-1)*Q*shock_vector;
end

%plot recent events and get share in GDP
figure
gov_spend_to_gdp_share=nom_gov_cons_inv./(real_gdp.*gdp_deflator/100);
plot(timeline,gov_spend_to_gdp_share)
axis tight
gov_spend_to_gdp_share=mean(gov_spend_to_gdp_share(timeline>1953.75 & timeline <2008));

figure
priv_cons_to_gdp_share=(nom_priv_cons_nd + nom_priv_cons_serv)./(real_gdp.*gdp_deflator/100);
plot(timeline,priv_cons_to_gdp_share)
axis tight
priv_cons_to_gdp_share=mean(priv_cons_to_gdp_share(timeline>1953.75 & timeline <2008));

%plot IRFs
figure
subplot(1,3,1);
plot(IRFs(G_pos,:));
title(var_names(G_pos));
ylim([-1 2.5])
subplot(1,3,2);
plot(1/gov_spend_to_gdp_share*IRFs(Y_pos,:));
title(var_names(Y_pos));
ylim([-1 2.5])

subplot(1,3,3);
plot(priv_cons_to_gdp_share/gov_spend_to_gdp_share*IRFs(C_pos,:));
title(var_names(C_pos));
ylim([-1 2.5])
