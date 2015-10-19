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

function [R] = cholpsd(X)
%%
% Upper Cholesky decomposition of psd matrix
% Tries to approach numerically $\lim_{\epsilon \to 0} chol(X+\epsilon I)$
%% Syntax
%   [R] = cholpsd(X)
%% Arguments
% * _X_ psd matrix _(n, n)_
%% Outputs
% * _R_ upper triangular matrix _(n, n)_ such that _X=R'R_
%% See also
% <solve_chol.html solve_chol>

m = min(min(X));
if m<0; fprintf('Error, X is negative in cholpsd\n'); error(''); end
m = max(eps, m*1e-14);
e = m;
I = eye(size(X));
ok = false;
while ~ok
    try
        R = chol(X);
        ok = true;
    catch
        % if the Cholesky decomposition failed, try to add a small epsilon on the diagonal
        X = X+e*I;
        if e > 1e6 * m
            fprintf('Warning, adding %f for cholpsd\n', e);
        end
        e = 10*e;
    end
end

end