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

function [nll] = gp_lik(Ktt, Yt, BayesInv, Ht)
%%
% Negative log likelihood given oberservations _Yt_ at _Xt_
%% Syntax
%   nll = gp_lik(Ktt, Yt, BayesInv)
%   nll = gp_lik(Ktt, Yt, BayesInv, Ht)
%% Arguments
% * _Ktt_ matrix _(nt, nt)_ of kernel between the points of _Xt_
% * _Yt_ vector _(nt, 1)_ of observations
% * _BayesInv_ structure array returned by _<gp_inf.html gp_inf>(Ht, Ktt, Yt, noise)_
% * _Ht_ matrix _(nt, b)_ of basis data as returned by _<basis_cst.html basis_cst>(Xt)_
%% Outputs
% * _nll_ float, negative of the logarithm of the likelihood
%% See also
% <gp_loolik.html gp_loolik> | <bfgs_search_prior.html bfgs_search_prior>

if nargin<4; Ht = []; end

n = size(Yt,1);
m = rank(Ht);
if m
    RK = cholpsd(Ktt);
    A = Ht'*(solve_chol(RK, Ht));
    RA = cholpsd(A);
    KHAHK = solve_chol(RK, Ht*solve_chol(RA, Ht'*solve_chol(RK, eye(n))));

    nll = .5*Yt'*BayesInv.invCY - .5*Yt'*KHAHK*Yt ...
          + sum(log(diag(RK))) + sum(log(diag(A))) + .5*(n-m)*log(2*pi);
else
    nll = .5*Yt'*(BayesInv.invCY) + .5*sum(log(diag(BayesInv.RC))) + .5*n*log(2*pi);
end

end
