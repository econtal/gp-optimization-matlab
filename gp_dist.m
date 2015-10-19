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

function [D2, Kuv_post] = gp_dist(Kuv, Ktu, Ktv, dKuu, dKvv, BayesInv, Ht, Hu, Hv)
%%
% Canonical GP distance $d^2(u,v) = Var[f(u)-f(v) \mid Xt] = \sigma_t^2(u)+\sigma_t^2(v)-2k(u,v)$
%% Syntax
%   D2 = gp_dist(Kuv, Ktu, Ktv, dKuu, dKvv, BayesInv)
%   D2 = gp_dist(Kuv, Ktu, Ktv, dKuu, dKvv, BayesInv, Ht, Hu, Hv)
%% Arguments
% * _Kuv_ kernel matrix _(nu, nv)_ between the points of _U_ and _V_
% * _Ktu_ kernel matrix _(nt, nu)_ between the points of _Xt_ and _U_
% * _Ktv_ kernel matrix _(nt, nv)_ between the points of _Xt_ and _V_
% * _dKuu_ vector _(nu, 1)_ of the diagonal kernel between the points of _U_
% * _dKvv_ vector _(nv, 1)_ of the diagonal kernel between the points of _V_
% * _BayesInv_ struct array as returned by _<gp_inf.html gp_inf>(Ht, Ktt, Yt, noise)_
% * _Ht_ matrix _(nt, b)_ basis for the points of _Xt_
% * _Hu_ matrix _(nu, b)_ basis for the points of _U_
% * _Hv_ matrix _(nv, b)_ basis for the points of _V_
%% Outputs
% * _D2_ matrix _(nu, nv)_ of squared distance between _U_ and _V_
%% See also
% <gp_pred.html gp_pred>

if nargin<7; Ht = []; end
if nargin<8; Hu = []; end
if nargin<9; Hv = []; end

[~,s2u] = gp_pred(Ktu, dKuu, BayesInv, Ht, Hu);
[~,s2v] = gp_pred(Ktv, dKvv, BayesInv, Ht, Hv);

% TODO update with basis
if Ht
    error('gp_dist is not implemented with basis functions.')
end
Kuv_post = Kuv - Ktu'*solve_chol(BayesInv.RC, Ktv);
D2 = bsxfun(@plus,s2u,bsxfun(@minus,s2v',2*Kuv_post));
D2(D2<eps) = 0;

end