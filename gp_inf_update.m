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

function BayesInv = gp_inf_update(Ktt12, Ktt22, Yt, noise, BayesInv1, Ht)
%%
% Bayesian system update with new observations at _Xt2_
%% Syntax
%   BayesInv = gp_inf_update(Ktt12, Ktt22, Yt, noise, BayInv1)
%   BayesInv = gp_inf_update(Ktt12, Ktt22, Yt, noise, BayInv1, Ht)
%% Arguments
% * _Ktt12_ matrix _(nt1, nt2)_ of kernel between the points of _Xt1_ and _Xt2_
% * _Ktt22_ matrix _(nt2, nt2)_ of kernel between the points of _Xt2_
% * _Yt_ vector _(nt1+nt2, 1)_ of observations
% * _noise_ noise standard deviation $\eta^2$
% * _BayInv1_ old system resolution for _Xt1_
% * _Ht_ matrix _(nt1+nt2, b)_ of basis data
%% Outputs
% struct array containing:
%
% * _RC_ upper triangular matrix _(nt1+nt2,nt1+nt2)_ of Cholesky decomposition of _Ktt+noise*I_
% * _invCY_ vector _(nt1+nt2,1)_ solution of $(K + \eta^2 I)^{-1} Y$
% * _beta_ vector _(b,1)_ solution of the basis system
%% See also
% <gp_inf.html gp_inf>

[nt1,nt2] = size(Ktt12);
if size(Ktt22,1)~=nt2 | size(Ktt22,2)~=nt2 | size(Yt,1)~=nt1+nt2
    error('Wrong arguments size: Ktt12:(%d,%d), Ktt22:(%d,%d), Yt:(%d,%d)',size(Ktt12),size(Ktt22),size(Yt));
end
if size(BayesInv1.RC,1)~=nt1
    error('Arguments size incompatible with previous BayesInv: Ktt12:(%d,%d), RC1:(%d,%d)',size(Ktt12),size(BayesInv1.RC));
end    

C12 = Ktt12;
if length(noise)==1
    C22 = Ktt22 + noise*eye(size(Ktt22)); % (nt x nt)
else
    C22 = Ktt22 + diag(noise); % (nt x nt)
end

% Cholsky updates (cf Osborne2010 p214)
RC1 = BayesInv1.RC;
RC12 = RC1'\C12;
RC22 = chol(C22 - RC12'*RC12);
RC = [RC1 RC12; zeros(size(RC12')) RC22]; % (nt x nt)

if Ht
    RHCH = cholpsd(Ht'*(solve_chol(RC, Ht))); % (b x b)

    % system resolution (cf RasmussenWilliams2006 Ch2 p28 Eq2.42)
    bet = solve_chol(RHCH, (Ht' * solve_chol(RC, Yt))); % (b x 1)
    invCY = solve_chol(RC, (Yt-Ht*bet));
else
    invCY = solve_chol(RC, Yt);
    bet = [];
    RHCH = [];
end

BayesInv = struct('RC',RC, 'invCY',invCY, 'bet',bet, 'RHCH',RHCH);

end
