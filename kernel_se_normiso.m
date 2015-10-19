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

function [K] = kernel_se_normiso(X, Y)
%%
% Isotropic squared exponential covariance function
%% Syntax
%   K = kernel_se_normiso(X,'diag')
%   K = kernel_se_normiso(X,Y)
%% Arguments
% * _X_ matrix _(nx, d)_ where _nx_ is the number of data points and _d_ is the dimension
% * _Y_ matrix _(ny, d)_ or 'diag' for diagonal self covariance 
%% Outputs
% * _K_ matrix _(nx, ny)_ or diagonal _(nx, 1)_
%% See also
% <kernel_se.html kernel_se>

ells = ones(size(X,2),1);
s = 1;

K = kernel_se(X, Y, s, ells);

end