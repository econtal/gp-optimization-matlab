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

function [f, Xs, Fs, Xt, Yt, Kss] = gp_sample(varargin)
%%
% Sample a Gaussian process in a hypercube and perform Bayesian inferance
%% Syntax
%   [f, Xs, Fs, Xt, Yt, Kss] = gp_sample(..., 'Name',Value)
%% Name-Value Pair Arguments
% * _d_ integer for the dimension (default: 1)
% * _size_ scalar for the length of the hypercube (default: 40)
% * _ns_ integer for the number of sampled points (default: 1000)
% * _nt_ integer for the number of training points (default: 10)
% * _kernel_ kernel function (default: <kernel_se_normiso.html @kernel_se_normiso>)
% * _basis_ basis function (default: <basis_none.html @basis_none>)
% * _noise_ scalar for the noise standard deviation (default: 0.01)
% * _plot_ boolean to visualize the optimization (default: true)
% * _posterior_ boolean to compute the Bayesian inferance (default: true)
% * _verbose_ boolean to monitor the process (default: true)
%% Outputs
% * _f_ function which returns noisy observations
% * _Xs_ matrix _(ns, d)_ of sampled data
% * _Fs_ matrix _(ns, 1)_ of sampled values
% * _Xt_ matrix _(nt, d)_ of training data
% * _Yt_ vector _(nt, 1)_ of training observations
% * _Kss_ matrix _(ns, ns)_ of kernel values


ip = inputParser;
ip.addOptional('d', 1);
ip.addOptional('size', 40);
ip.addOptional('ns', 1000);
ip.addOptional('nt', 10);
ip.addOptional('kernel', @kernel_se_normiso);
ip.addOptional('basis', @basis_none);
ip.addOptional('noise', 1e-2);
ip.addOptional('plot', true);
ip.addOptional('posterior', true);
ip.addOptional('verbose', true);
ip.parse(varargin{:});
opt = ip.Results;

if opt.verbose; fprintf('generating fn...\n'); end
Xs = opt.size * rand(opt.ns, opt.d);
Xt = Xs(1:opt.nt,:);
Kss = opt.kernel(Xs,Xs);
Fs = cholpsd(Kss)' * randn(opt.ns, 1);
f = @(X) Fs(X) + opt.noise * randn(size(X,1), 1);
Yt = f(1:opt.nt);

if opt.posterior
    if opt.plot; plot(Xt, Yt, 'r+'); hold on; end

    if opt.verbose; fprintf('Bayesian inferance...\n'); end
    Ktt = Kss(1:opt.nt,1:opt.nt);
    Ht = opt.basis(Xt);
    BayesInv = gp_inf(Ktt, Yt, opt.noise, Ht);
    
    if opt.verbose; fprintf('Prediction...\n'); end
    Kts = Kss(1:opt.nt,:);
    Hs = opt.basis(Xs);
    dKss = diag(Kss);
    [mu, s2] = gp_pred(Kts, dKss, BayesInv, Ht, Hs);

    if opt.plot
        [Z,IZ] = sort(Xs(:,1));
        plot(Z, mu(IZ), 'k');
        fill([Z; flipdim(Z,1)], [mu(IZ)+2*s2(IZ); flipdim(mu(IZ)-2*s2(IZ),1)], ...
             [7 8 7]/8, 'EdgeColor','None', 'FaceAlpha',.8);
    end
else
    if opt.plot
        if opt.d == 1
            [Z,IZ] = sort(Xs(:,1));
            plot(Z, Fs(IZ), 'k');
        elseif opt.d == 2
            scatter(Xs(:,1),Xs(:,2),30,Fs,'fill');
        end
    end
end

end