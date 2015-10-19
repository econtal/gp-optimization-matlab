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

function [eNet, PiX, DeltaT, U] = chaining_tree(D2, dmax, step, iMin, iMax, varargin)
%%
% Compute the chaining tree given the canonical distance matrix between point of _X_
%% Syntax
%   [eNet, PiX, DeltaT, U] = chaining_tree(D2, dmax, step, iMin, iMax)
%   [eNet, PiX, DeltaT, U] = chaining_tree(..., 'Name',Value)
%% Arguments
% * _D2_ matrix _(n, n)_ of canonical squared distance
% * _dmax_ scalar of initial diameter of _X_
% * _step_ scalar > 1 for the geometric decay of $\epsilon_i$
% * _iMin_ integer for the first level to consider
% * _iMax_ integer for the last level to consider
%% Name-Value Pair Arguments
% * _a_ scalar > 1 power to use for the geometric decay in the union bound, e.g. 2
% * _lza_ scalar of logarithm of the Riemann zeta of a, e.g. _log(pi^2/6)_
%% Outputs
% * _eNet_ vector _(1, N)_ of indices of points in the final $\epsilon$-net
% * _PiX_ matrix _(h, n)_ of indices of the closest element of the net for all _h=iMax-iMin_ levels
% * _DeltaT_ matrix _(h, N)_ of diameters of the cells of the net for all _h_ levels
% * _U_ vector _(h, 1)_ of negative log probabilities w.r.t the union bounds for all _h_ levels
%% See also
% <gp_dist.html gp_dist> | <enet_greedy.html enet_greedy>

ip = inputParser;
ip.addOptional('a', 2.2);
ip.addOptional('lza', 0.399141);
ip.parse(varargin{:});
opt = ip.Results;

n = size(D2,1);

Ti = [];
nSteps = iMax-iMin;
PiX = zeros(nSteps,n,'uint32');
DeltaT = inf(nSteps,n);
U = inf(nSteps,1);

for ind=1:nSteps
    i = ind+iMin-1;
    ei = dmax*step^-i;
    [ni,Ti] = enet_greedy(D2,ei,Ti);
    U(ind) = log(ni+1) + opt.a*log(i) + opt.lza;

    [TiDist2, TiInd] = min(D2(Ti,:), [], 1);
    PiX(ind,:) = TiInd;
    for j=1:length(Ti)
        DeltaTi_j = sqrt(max(TiDist2(TiInd==j)));
        if ~isempty(DeltaTi_j)
            DeltaT(ind,j) = DeltaTi_j;
        end
    end
end

eNet = Ti;

end