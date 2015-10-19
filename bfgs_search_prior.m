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

function [HP, kernel, noise, fval] = bfgs_search_prior(Xt, Yt, HPini, kfun, basis)
%%
% Optimize prior hyper-parameter with respect to the pseudo-likelihood, using the BFGS algorithm
%% Syntax
%   HP = bfgs_search_prior(Xt, Yt, HPini, kfun, basis)
%% Arguments
% * _Xt_ matrix _(n, d)_ where _n_ is the number of data points and _d_ is the dimension
% * _Yt_ vector _(n, 1)_ of noisy observations
% * _HPini_ vector _(h, 1)_ of initial hyper-parameters formatted as : _[log(sf2) log(sn) log(w1) log(w2)]_
% * _kfun_ kernel function such as _<kernel_se.html kernel_se>_
% * _basis_ basis function such as _<basis_none.html basis_none>_
%% Outputs
% * _HP_ vector _(h, 1)_ of locally optimal log hyper-parameters
% * _kernel_ found kernel function
% * _noise_ found noise standard deviation
%% See also
% <gp_loolik.html gp_loolik>

NelderMeadIters = 50;

Ht = basis(Xt);
f = @(hp) nllcost(Xt, Yt, hp(1), hp(2), hp(3:end), kfun, Ht);

% Starts with few Nelder-Mead iterations
[hpopt, fvalNM] = fminsearch(f, HPini, optimset('maxFunEvals',NelderMeadIters,'display','off'));
% BGFS search
[hpopt, fval] = fminunc(f, hpopt,optimset('display','off','LargeScale','off'));

HP = hpopt;
kernel = @(x,y) kfun(x,y, exp(HP(1)), exp(HP(3:end)));
noise = exp(HP(1)+HP(2));

end


function [nll] = nllcost(Xt, Yt, sf2, rsn, W, kfun, Ht)

Ktt = kfun(Xt, Xt, exp(sf2), exp(W));
noise = exp(sf2+rsn);
BayesInv = gp_inf(Ktt, Yt, noise, Ht);
nll = gp_loolik(Ktt, Yt, noise, BayesInv, Ht);

end