function [minus_log_posterior,alarm_flag,F,H,R,Q]=posterior_Kalman(par_vector,y)
%function [minus_log_posterior,alarm_flag,F,H,R,Q]=posterior_Kalman(par_vector,y)
% computes the log-likelihood of a DSGE model using the Kalman filter
% Inputs:
%   par_vector:     [npar_estimate by 1] vector of parameter over which to optimize
%   Y:              [nobservables by T] Matrix of observed variables
%
% Output:
%   minus_log_posterior [scalar] value of minus the log-posterior function
%   alarm_flag:         [scalar] 1 if model could not be solved, 0 otherwise
%   F                   [nstates by nstates] state transition matrix
%   Q                   [nstates by nshocks] covariance matrix of structural
%                                           shocks
%   H                   [n_observables by nstates] observation to states mapping matrix
%   R                   [n_observables by n_measurement_errors] covariance
%                                       matrix of measurement errors

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

alarm_flag=0;

if size(y,1)>size(y,2)
    y=y';
end
[nobservables, T]=size(y);
% rebuild parameter vector
[rhoz, rhog, eps_z, eps_g]=par_retransform_bounded(par_vector(1),par_vector(2),par_vector(3),par_vector(4));

[priorval,alarm_flag_1]=evaluate_prior([rhoz, rhog, eps_z, eps_g]);
if alarm_flag_1
    minus_log_posterior=Inf;
    alarm_flag=2;
else
    
    [F,H,alarm_flag_2]=get_solution_RBC(rhoz,rhog);
    H=[H;
        0 0 1];
    if alarm_flag_2
        minus_log_posterior=10^8+sum(par_vector.^2); %penalty function to avoid discontinuities
        alarm_flag=1;
    else
        %order x=[k_t z_t g_t]; y=[c_t]
        nstates=length(F);
        R=zeros(nobservables); %no measurement error
        %build AR matrices
        Q_chol=diag([0 eps_z, eps_g]); %set variance of shock to k_t to 0 and the others to estimated parameters
        Q=Q_chol*Q_chol';
        
        %Matrices in Kalman filter (follows Hamilton 1994):
        % - Sigmatgtm1: MSE of x_t given y_{t-1}
        % - xhattgtm1: Forecast of x_{t} given y_{t-1}
        % - Sigmatp1gt: MSE of x_{t+1} given y_{t}
        % - xhattp1gt: Forecast of x_{t+1} given y_{t}
        % - K: Kalman gain
        
        %initialize variables at unconditional mean and variance
        Sigmatgtm1 = lyapunov_symm(F,Q,1e-10,1e-6,1e-15);
        xhattgtm1=zeros(length(F),1);
        %% Kalman Filter Loop
        loglik = 0; % initialize log-likelihood
        %% Kalman Filter Loop
        for t=1:T
            % gain matrix
            omega=H*Sigmatgtm1*H'+R; %13.2.25
            if det(omega)<=0
                minus_log_posterior=10^8+sum(par_vector.^2); %penalty function to avoid discontinuities
                return
            end
            K=F*Sigmatgtm1*H'/omega; %13.2.19
            %Calculate Zeta hat and P at t and t+1 given info a t
            dataresid=y(:,t)-H*xhattgtm1; %part of 13.2.20
            xhattp1gt=F*xhattgtm1+K*dataresid; %13.2.20
            Sigmatp1gt=(F-K*H)*Sigmatgtm1*(F'-H'*K')+K*R*K'+Q; %13.2.28
            %want to minimize, take opposite sign
            logliktemp=-nobservables/2*log(2*pi)-0.5*reallog(det(omega))+((-0.5*dataresid'/omega*dataresid));
            loglik=logliktemp+loglik;
            %reset some values for next step
            Sigmatgtm1=Sigmatp1gt;
            xhattgtm1=xhattp1gt;
        end
        minus_log_posterior =-(loglik+priorval);
    end
end