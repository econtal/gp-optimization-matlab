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

function [BayesInv] = gp_inf(Ktt, Yt, noise, Ht)
%%
% Bayesian system resolution for computing posterior of GP given the observations _Yt_ at _Xt_
%% Syntax
%   BayesInv = gp_inf(Ktt, Yt, noise)
%   BayesInv = gp_inf(Ktt, Yt, noise, Ht)
%% Arguments
% * _Ktt_ matrix _(nt, nt)_ of kernel between the points of _Xt_
% * _Yt_ vector _(nt, 1)_ of observations
% * _noise_ noise standard deviation $\eta^2$
% * _Ht_ matrix _(nt, b)_ of basis data as returned by _<basis_cst.html basis_cst>(Xt)_
%% Outputs
% struct array containing:
%
% * _RC_ upper triangular matrix _(nt,nt)_ of Cholesky decomposition of _Ktt+noise*I_
% * _invCY_ vector _(nt,1)_ solution of $(K + \eta^2 I)^{-1} Y$
% * _beta_ vector _(b,1)_ solution of the basis system
%% See also
% <gp_pred.html gp_pred> | <gp_inf_update.html gp_inf_update>

if nargin<4; Ht = []; end

C = Ktt + noise*eye(size(Ktt)); % (nt x nt)

% Cholesky decomposition
RC = cholpsd(C); % (nt x nt)
if Ht
    RHCH = cholpsd(Ht' * (solve_chol(RC, Ht))); % (b x b)

    % system resolution
    bet = solve_chol(RHCH, (Ht' * solve_chol(RC, Yt))); % (b x 1)
    invCY = solve_chol(RC, (Yt - Ht*bet));
else
    invCY = solve_chol(RC, Yt);
    bet = [];
    RHCH = [];
end

BayesInv = struct('RC',RC, 'invCY',invCY, 'bet',bet, 'RHCH',RHCH);

end
