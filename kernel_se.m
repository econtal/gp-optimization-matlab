% This file is part of GpOptimization.
%
% GpOptimization is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% GpOptimization is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with GpOptimization. If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (c) by Emile Contal, 2015

function [K] = kernel_se(X, Y, s, ells)
%%
% Squared exponential covariance function
%% Syntax
%   K = kernel_se(X,'diag',s,ells)
%   K = kernel_se(X,Y,s,ells)
%% Arguments
% * _X_ matrix _(nx, d)_ where _nx_ is the number of data points and _d_ is the dimension
% * _Y_ matrix _(ny, d)_ or 'diag' for diagonal self covariance 
% * _s_ scalar _(1,1)_ for covariance scale
% * _ell_ vector _(1,d)_ for covariance length-scales
%% Outputs
% * _K_ matrix _(nx, ny)_ or diagonal _(nx, 1)_
%% See also
% <kernel_matern.html kernel_matern> | <kernel_se.html kernel_se>

ARD = diag(1./ells);
if isa(Y, 'char') && strcmp(Y, 'diag')
    D = zeros(size(X,1),1);
else
    D = sq_dist(ARD*X', ARD*Y');
end

K = s * exp(-D/2);

end