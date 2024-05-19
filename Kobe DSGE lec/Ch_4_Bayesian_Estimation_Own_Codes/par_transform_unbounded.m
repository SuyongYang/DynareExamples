function [rhoz, rhog, eps_z, eps_g]=par_transform_unbounded(rhoz,rhog,eps_z,eps_g)
% [rhoz, rhog, eps_z, eps_g]=par_transform_unbounded(rhoz,rhog,eps_z,eps_g)
% transforms parameters with bounded support into parameters with unbounded
% support by using either exp() or logistic transformation;

% for parameters with support [LB Inf): log(x-LB)
% for parameters with support [LB UB]: log((UB-x)./(x-LB));

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

eps_z=log(eps_z);
eps_g=log(eps_g);
rhoz= log((0.9999-rhoz)./(rhoz-(-0.9999)));
rhog= log((0.9999-rhog)./(rhog-(-0.9999)));

