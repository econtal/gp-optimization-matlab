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

function M = cummax(R)
%%
% Cumulative maximum as used in simple regrets
%% Syntax
%   M = cummax(R)
%% Arguments
% * _R_ matrix _(n,t)_ for _n_ run of length _t_
%% Outputs
% * _M_ matrix _(n,t)_ of maximum along the column _(:,1:i)_ for each row

[d,n] = size(R);
M = zeros(d,n);
m = R(:,1);
M(:,1) = m;

for i=2:n
    m = max(m, R(:,i));
    M(:,i) = m;
end

end