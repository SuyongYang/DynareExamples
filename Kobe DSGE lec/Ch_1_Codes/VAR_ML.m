function [minus_loglik,Sigma_U]=VAR_ML(par_vector,Y,Z)
%function [minus_loglik,Sigma_U]=VAR_ML(par_vector,Y,Z)
% computes minus the log-likelihood of a VAR
% Inputs:
%   par_vector:     [(constant_dummy+trend_dummy)*nvars+nvars^2*nlags+nvars*(nvars-1)/2 by 1] vector of parameter over which to optimize
%   Y:              [nvars by T] Matrix of observed variables 
%   Z:              [nvars*nlags by T] Matrix of regressors
%
% Outputs:
%   minus_loglik    [scalar]        negative of the log-likelihood function
%   Sigma_U         [nvars*nvars]   covariance matrix of the shocks

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

if size(Y,1)>size(Y,2)
   error('Wrong dimension of Y') 
end
[nvars, T]=size(Y);
Y=Y(:);

%build AR matrices
betta=par_vector(1:end-nvars*(nvars+1)/2,1);
L=vec2cholmat(par_vector(end-nvars*(nvars+1)/2+1:end,1));

Sigma_U=L*L';
loglik=-nvars*T/2*log(2*pi)-T/2*log(det(Sigma_U))...
    -0.5*(Y-kron(Z',eye(nvars))*betta)'*kron(eye(T),eye(nvars)/Sigma_U)*...
    (Y-kron(Z',eye(nvars))*betta);
minus_loglik=-loglik;