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

function [mu, sigma2] = gp_downdate(Ktt, Yt, i, BayesInv, Ht)
%%
% Posterior $\mu(x_i)$ and $\sigma^2(x_i)$ given $X_t \setminus \{x_i\}$
%% Syntax
%   [mu, sigma2] = gp_downdate(Ktt, Yt, i, BayesInv)
%   [mu, sigma2] = gp_downdate(Ktt, Yt, i, BayesInv, Ht)
%% Arguments
% * _Ktt_ kernel matrix _(nt, nt)_ between the points of _Xt_
% * _Yt_ vector _(nt, 1)_ of observations
% * _i_ indice of removed observation
% * _BayesInv_ struct array returned by _<gp_inf.html gp_inf>(Ht, Ktt, Yt, noise)_
% * _Ht_ matrix _(nt, b)_ of basis data as returned by _<basis_cst.html basis_cst>(Xt)_
%% Outputs
% * _mu_ scalar _(1,1)_ posterior mean $E[f(x_i) \mid X_t\setminus x_i, Y_t \setminus y_i]$
% * _sigma2_ scalar _(1,1)_ posterior variance $V[f(x_i) \mid X_t\setminus x_i, Y_t \setminus y_i]$
%% See also
% <gp_pred.html gp_pred> | <gp_loolik.html gp_loolik>

n = size(Ktt,1);
T = [1:i-1 i+1:n];
Yt1 = Yt(T);

% Covariance
Kti = Ktt(T,i);
Kii = Ktt(i,i);

% Cholsky downdates (cf Osborne2010 p216)
RC = BayesInv.RC;
RC11 = RC(1:i-1,1:i-1);
RC13 = RC(1:i-1,i+1:end);
S23  = RC(i,i+1:end);
S33  = RC(i+1:end,i+1:end);

RC33 = cholupdate(S33, S23');
RC1 = [RC11 RC13; zeros(size(RC13')) RC33];

if Ht
    Ht1 = Ht(T,:);
    Hi = Ht(i,:);

    RHCH1 = cholpsd(Ht1'*solve_chol(RC1,Ht1));

    % System resolution (cf RasmussenWilliams2006 Ch2 p28 Eq2.42)
    Ri = Hi - Ht1'*solve_chol(RC1,Kti);
    bet = solve_chol(RHCH1, (Ht1'*solve_chol(RC1, Yt1)));
    invCY = solve_chol(RC1, (Yt1 - Ht1*bet));
    mu = Hi'*bet + Kti'*invCY;
else
    invCY = solve_chol(RC1, Yt1);
    mu = Kti'*invCY;
    bet = [];
end

% sigma2
if nargout > 1
    Vf = RC1'\Kti;
    covf = Kii - sum(Vf.*Vf, 1)';

    if Ht
        Vb = RHCH1'\Ri;
        covb = sum(Vb.*Vb, 1)';
        sigma2 = covb + covf;
    else
        sigma2 = covf;
    end
end

end
