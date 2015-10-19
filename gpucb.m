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

function [ucb] = gpucb(s2, u, n)
%%
% Compute the ucb as in the GP-UCB algorithm
%% Syntax
%   ucb = chaining_ucb(s2, u, n)
%% Arguments
% * _s2_ matrix _(n, 1)_ of posterior variance
% * _u_ scalar for negative log probability
% * _n_ number of test points
%% Outputs
% * _ucb_ vector _(n, 1)_ such that $P[\forall x,~f(x)-\mu(x)>ucb(x)] < \exp(-u)$
%% See also
% <chaining_ucb.html chaining_ucb>

ucb = sqrt(2*(u+log(n))*s2);

end