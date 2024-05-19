function [log_priorval,alarm]=evaluate_prior(par)
% [priorval,alarm]=evaluate_prior(par)
% evaluates the prior PDF at the parameters given in par
% Inputs:
%   - par   [npar*1]    vector of parameters
% 
% Ouputs:
%   - log_priorval  [scalar]     log prior density
%   - alarm         [scalar]     indicator that is 1 if there is a violation of the prior bounds

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


alarm=0;
log_priorval=0; % intialize prior

% rhoz,rhog
[a, b]=betaparametertransform(0.7,0.1^2);
ii=[1:2];
prioradd=log(betapdf(0.999^(-1)*par(ii)',a,b)); % 1/0.999 to ensure non-unit root
if ~isinf(prioradd)
  log_priorval=log_priorval+sum(prioradd);
else
  alarm=1;  
  return
end

% eps_z, eps_g
[a,b] = inverse_gamma_specification(0.01,0.1,1,0);
ii=[3:4];
prioradd=sum(lpdfig1(par(ii)',a,b));
if ~isinf(prioradd)
  log_priorval=log_priorval+prioradd;
else
  alarm=1;  
  return
end

end
