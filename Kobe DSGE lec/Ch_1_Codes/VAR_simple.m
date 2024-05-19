% Copyright (C) 2014-2016 Benjamin Born and Johannes Pfeifer
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

%% simulate data

n_obs = 10000;
x = zeros(2,n_obs);
cov_mat = [1 0.02; 0.02 1]; % eye(2)

chol_cov = chol(cov_mat,'lower');
A_mat = [0.5 0.02; 0.01 0.7];

epsilon = randn(2,n_obs);

for ii = 2:n_obs
    x(:,ii) = A_mat*x(:,ii-1) + chol_cov*epsilon(:,ii);
end

%% Estimate VAR

Ymat = x(:,2:end)';
Zmat = x(:,1:end-1)';

[nobs,nvars] = size(Ymat);

Z = Zmat';
ymat = Ymat';

bhat = kron((Z*Z')\Z,eye(nvars))*ymat(:);