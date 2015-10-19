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

function [N, eNet] = enet_greedy (D2, e, eNetInit)
%%
% Compute an $\epsilon$-net of _Xt_ given the distance matrix
%% Syntax
%   [N, eNet] = enet_greedy(D2, e, eNetInit)
%% Arguments
% * _D2_ matrix _(n, n)_ of canonical squared distance
% * _e_ scalar for $\epsilon$
% * _eNetInit_ vector _(1, N)_ of indices of points in the initial $\epsilon$-net
%% Outputs
% * _N_ scalar for the size of the final $\epsilon$-net
% * _eNet_ vector _(1, N)_ of indices of points in the final $\epsilon$-net
%% See also
% <gp_dist.html gp_dist> | <chaining_tree.html chaining_tree>

e2 = e^2;
[n,~] = size(D2);
L = sparse(D2<=e2);
eNet = eNetInit;
if eNet
    nc = ~any(L(:,eNet),2); % not covered
else
    nc = true(n,1);
end
while any(nc)
    nc_ind = find(nc)';
    P = sum(L(nc,nc), 2); % number of points in the ball of center X_i and radius e
    [~,maxP] = max(P);                    % we choose the maximum
    isolated = find(P==1)';               % and all the isolated points
    eNet = [eNet nc_ind([maxP isolated])];
    nc(nc) = ~any(L(nc,eNet),2);
end
N = length(eNet);

end
