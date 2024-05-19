function [a,b]=betaparametertransform(betamean,betavariance)
% Compute the location and shape parameter of the Beta distribution from the mean and variance
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

b=(betamean-2*betamean^2+betamean^3-betavariance+betamean*betavariance)/betavariance;
a=-betamean*b/(betamean-1);