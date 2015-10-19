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

function [ucb] = chaining_ucb(D2, s2, u, varargin)
%%
% Compute the chaining UCB given the canonical distance matrix between point of _X_
%% Syntax
%   ucb = chaining_ucb(D2, s2, u)
%   ucb = chaining_ucb(..., 'Name',Value)
%% Arguments
% * _D2_ matrix _(n, n)_ of canonical squared distance
% * _s2_ matrix _(n, 1)_ of posterior variance
% * _u_ scalar for negative log probability
%% Name-Value Pair Arguments
% * _step_ scalar > 1 power to use for the geometric decay of $\epsilon$
% * _cdfinv_ function to upper bound the tail probabilities e.g. _@(u,d) d.*sqrt(2*u)_
%% Ouputs
% * _ucb_ vector _(n, 1)_ such that $P[\sup_{x^\star}f(x^\star)-f(x)-\mu(x^\star)+\mu(x^\star) > ucb(x)+cst] < \exp(-u)$
%% See also
% <gpucb.html gpucb> | <chaining_tree.html chaining_tree>

ip = inputParser;
ip.addOptional('step', 2);
ip.addOptional('cdfinv', @(u,d) sqrt(2)*d.^2*erfcinv(2*exp(-u)));
ip.parse(varargin{:});
opt = ip.Results;

n = size(D2,1);

% chaining tree
dmax = max(sqrt(2*s2));
deltaInd = ceil(-log(sqrt(s2)./dmax)/log(opt.step)); % indices such that an upper bound on Delta(X) is < s(X)
[Tree,PiX,DeltaT,U] = chaining_tree(D2,dmax,opt.step,min(deltaInd),max(deltaInd));

% compute ucb(X)
ucb = zeros(n,1);
L1 = false(n,1); % i > min{ i: Delta_i(X) < s(X) }
DeltaXi1 = inf(n,1);
for i=1:size(PiX,1)
    ui = U(i)+u;
    DeltaXi = DeltaT(i,PiX(i,:))';
    LT = DeltaXi < sqrt(s2);
    L0 = LT & (~L1);
    ucb(L0) = opt.cdfinv(ui, sqrt(s2(Tree(PiX(i,L0))))); % i = min{ i: Delta_i(X) < s(X) }
    ucb(L1) = ucb(L1) + opt.cdfinv(ui, DeltaXi1(L1)); % i > min{ }
    DeltaXi1 = DeltaXi;
    L1(LT) = true;
end

end