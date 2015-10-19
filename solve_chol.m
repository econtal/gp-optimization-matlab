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

function [X] = solve_chol(R, Y)
%%
% Linear system resolution $MX=Y$ given upper Cholesky decomposition of $M=R^\top R$
%% Syntax
%   [X] = solve_chol(R, Y)
%% Arguments
% * _R_ upper triangular matrix _(n, n)_ of Cholesky decomposition of _M_
% * _Y_ matrix _(n, 1)_
%% Outputs
% * _X_ matrix _(n, 1)_ such that _M*X=Y_
%% See also
% <cholpsd.html cholpsd>

X = R\(R'\Y);