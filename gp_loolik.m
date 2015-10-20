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

function [nll] = gp_loolik(Ktt, Yt, BayesInv, Ht)
%%
% Negative log leave-one-out likelihood also called Pseudo-likelihood given oberservations _Yt_ at _Xt_
%% Syntax
%   nll = gp_loolik(Ktt, Yt, BayesInv)
%   nll = gp_loolik(Ktt, Yt, BayesInv, Ht)
%% Arguments
% * _Ktt_ matrix _(nt, nt)_ of kernel between the points of _Xt_
% * _Yt_ vector _(nt, 1)_ of observations
% * _BayesInv_ structure array returned by _<gp_inf.html gp_inf>(Ht, Ktt, Yt, noise)_
% * _Ht_ matrix _(nt, b)_ of basis data as returned by _<basis_cst.html basis_cst>(Xt)_
%% Output
% * _nll_ float, negative of the logarithm of the pseudo-likelihood
%% See also
% <gp_lik.html gp_lik> | <bfgs_search_prior.html bfgs_search_prior>

if nargin<5; Ht = []; end

n = size(Ktt, 1);
nll = 0;

for i=1:n
    [mui, s2i] = gp_downdate(Ktt, Yt, i, BayesInv, Ht);
    nll = nll + .5*log(s2i) + (Yt(i)-mui)^2/(2*s2i) + .5*log(2*pi);
end

if ~isreal(nll)
    nll = inf;
end

end
