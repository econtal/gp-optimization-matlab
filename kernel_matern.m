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

function [K] = kernel_matern(X, Y, s, ells, nu)
%%
% Matern covariance function for $\nu \in \{ \frac 1 2, \frac 3 2, \frac 5 2\}$
%% Syntax
%   K = kernel_matern(X,'diag',s,ells,nu)
%   K = kernel_matern(X,Y,s,ells,nu)
%% Arguments
% * _X_ matrix _(nx, d)_ where _nx_ is the number of data points and _d_ is the dimension
% * _Y_ matrix _(ny, d)_ or 'diag' for diagonal self covariance 
% * _s_ scalar _(1,1)_ for covariance scale
% * _ells_ vector _(1,d)_ of covariance length-scales
% * _nu_ scalar _(1,1)_ for Matern parameter $2\nu$, available values are _1_,_3_ and _5_
%% Outputs
% * _K_ matrix _(nx, ny)_ or diagonal _(nx, 1)_
%% See also
% <kernel_se.html kernel_se>

ARD = diag(1./ells);

if isa(Y, 'char') && strcmp(Y, 'diag')
    D = zeros(size(X,1),1);
else
    D = sqrt(sq_dist(sqrt(nu)*ARD*X', sqrt(nu)*ARD*Y'));
end

switch nu
  case 1, K = s * exp(-D);
  case 3, K = s * (1+D) .* exp(-D);
  case 5, K = s * (1+D.*(1+D/3)) .* exp(-D);
  otherwise
    error('kernel_matern is only defined for n=1,3 or 5')
end

end