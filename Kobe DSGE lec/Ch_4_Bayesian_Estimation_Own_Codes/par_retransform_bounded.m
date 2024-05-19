function  [rhoz, rhog, eps_z, eps_g]=par_retransform_bounded(rhoz,rhog,eps_z,eps_g);
% [rhoz, rhog, eps_z, eps_g]=par_retransform_bounded(rhoz,rhog,eps_z,eps_g);
% retransforms parameters with unbounded support to their origina values with bounded
% support 
% for parameters with support [LB Inf): LB+exp(x)   (inverse transformation of log(x-LB))
% for parameters with support [LB UB]: (UB+exp(x)*LB)./(1+exp(x)) (inverse tranformation of log((UB-x)./(x-LB)))

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

eps_z=exp(eps_z);
eps_g=exp(eps_g);
rhoz= (0.9999+exp(rhoz)*(-0.9999))./(1+exp(rhoz));
rhog= (0.9999+exp(rhog)*(-0.9999))./(1+exp(rhog));

