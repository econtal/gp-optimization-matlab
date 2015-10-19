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

function [queries, Yt] = gpopt(f, Xt, Yt, Xs, T, varargin)
%%
% Sequential optimization using GP algorithms
%% Syntax
%   [queries, Yt] = gpopt(f, Xt, Yt, Xs, T)
%   [queries, Yt] = gpopt(..., 'Name',Value)
%% Arguments
% * _f_ function which returns observations
% * _Xt_ matrix _(nt, d)_ of initial data
% * _Yt_ vector _(nt, 1)_ of initial observations
% * _Xs_ matrix _(ns, d)_ of search data
% * _T_ integer for the number of iterations
%% Name-Value Pair Arguments
% * _algo_ string code for the query algorithm
%    possible values are 'gpucb' (default), 'chaining', 'ei' for purely sequential optimization
%    and 'gpucbpe', 'greedyucb' for batch sequential optimization
% * _noise_ scalar for the noise standard deviation (default: 1e-2)
% * _u_ scalar for the negative logarithm of the upper bound probability (default: 3)
% * _B_ integer for the size of the batch in batch sequential optimization (default: 1)
% * _kernel_ kernel function for the inferance (default: <kernel_se_normiso.html @kernel_se_normiso>)
% * _basis_ basis function for the inferance (default: <basis_none.html @basis_none>)
% * _Kss_ matrix _(ns,ns)_ of the kernel between points of _Xs_, used by the 'chaining' algorithm
% * _plot_ boolean to visualize the optimization (default: false)
% * _verbose_ boolean to monitor the optimization (default: true)
%% Outputs
% * _queries_ vector _(1,T)_ of queries indices
% * _Yt_ vector _(T,1)_ of observations
%% See also
% <gpucb.html gpucb> | <chaining_ucb.html chaining_ucb>

ip = inputParser;
ip.addOptional('algo', 'gpucb');
ip.addOptional('u', 3);
ip.addOptional('B', 1);
ip.addOptional('kernel', @kernel_se_normiso);
ip.addOptional('basis', @basis_none);
ip.addOptional('noise', 1e-2);
ip.addOptional('plot', false);
ip.addOptional('verbose', true);
ip.addOptional('Kss', []);
ip.parse(varargin{:});
opt = ip.Results;

% check input
if strcmp(opt.algo, 'chaining') & isempty(opt.Kss)
    error('You must provide Kss in the optional arguments for using the chaining algorithm');
end
if opt.B>1 & ~(strcmp(opt.algo,'gpucbpe') | strcmp(opt.algo,'gpbucb') | strcmp(opt.algo,'greedyucb'))
    error('The algorithm is not adapted for batch optimization');
end
target_plot = [];

[ns,d] = size(Xs);

queries = [];
Ht = opt.basis(Xt);
Hs = opt.basis(Xs);
Ktt = opt.kernel(Xt,Xt);
Kts = opt.kernel(Xt,Xs);
dKss = opt.kernel(Xs,'diag');

BayesInv = gp_inf(Ktt, Yt, opt.noise, Ht);

for iter=1:T
    % GP pred
    [mu, s2] = gp_pred(Kts, dKss, BayesInv, Ht, Hs);

    % compute ucb
    u = opt.u + log(iter^2*pi^2/6);
    switch opt.algo
      case 'gpucb'
        % GP-UCB - Srinivas et al. (2012)
        ucb = gpucb(s2, u, ns);
        target = mu + ucb;
        [target_val,xt] = max(target);
      case 'chaining'
        % Chainig-UCB - Contal et al. (2015)
        D2 = gp_dist(opt.Kss, Kts, Kts, dKss, dKss, BayesInv);
        ucb = chaining_ucb(D2, s2, opt.u + log(iter^2*pi^2/6));
        target = mu + ucb;
        [target_val,xt] = max(target);
      case 'ei'
        % Expected Improvement
        fmax = max(Yt);
        ni = (mu-fmax) ./ sqrt(s2);
        target = (mu-fmax).*normcdf(ni) + sqrt(s2) .* normpdf(ni);
        [target_val,xt] = max(target);
      case 'greedyucb'
        % Batch-greedy ucb policy, similar to GP-BUCB without the initialization phase
        xt = [];
        for b=1:opt.B
            ucb = gpucb(s2, u, ns);
            target = mu + ucb;
            [target_val,xtb] = max(target);
            xt = [xt; xtb];
            if b<opt.B
                Ht = [Ht; opt.basis(Xs(xtb,:))];
                Ktt12 = opt.kernel(Xt,Xs(xtb,:));
                Ktt22 = opt.kernel(Xs(xtb,:),Xs(xtb,:));
                BayesInv = gp_inf_update(Ktt12, Ktt22, [Yt; mu(xt)], opt.noise, BayesInv, Ht);
                Kts = [Kts; opt.kernel(Xs(xtb,:),Xs)];
                [~, s2] = gp_pred(Kts, dKss, BayesInv, Ht, Hs);
                Xt = [Xt; Xs(xtb,:)];
            end
        end
      case 'gpucbpe'
        % GP-UCB-PE - Contal et al. (2013)
        ucb = gpucb(s2, u, ns);
        target = mu + ucb;
        [target_val,xtb] = max(target);
        xt = [xtb];
        lcb = mu - ucb;
        ydot = max(lcb);
        is_in_Rt = target >= ydot;
        for b=2:opt.B
            % update
            Ht = [Ht; opt.basis(Xs(xtb,:))];
            Ktt12 = opt.kernel(Xt,Xs(xtb,:));
            Ktt22 = opt.kernel(Xs(xtb,:),Xs(xtb,:));
            BayesInv = gp_inf_update(Ktt12, Ktt22, [Yt; mu(xt)], opt.noise, BayesInv, Ht);
            Kts = [Kts; opt.kernel(Xs(xtb,:),Xs)];
            [~, s2] = gp_pred(Kts, dKss, BayesInv, Ht, Hs);
            Xt = [Xt; Xs(xtb,:)];
            % next query
            s2_in_Rt = s2 .* is_in_Rt;
            [~, xtb] = max(s2_in_Rt);
            xt = [xt; xtb];
        end
    end
    
    % query point
    queries = [queries xt];
    yt = f(xt);
    Xt = [Xt; Xs(xt(end),:)]; % if batch, Xt is already partially updated
    Yt = [Yt; yt];

    % GP update
    Ht = [Ht; opt.basis(Xt(end,:))];
    Ktt12 = opt.kernel(Xt(1:end-1,:),Xt(end,:));
    Ktt22 = opt.kernel(Xt(end,:),Xt(end,:));
    BayesInv = gp_inf_update(Ktt12, Ktt22, Yt, opt.noise, BayesInv, Ht);
    Kts = [Kts; opt.kernel(Xt(end,:),Xs)];

    % monitoring
    switch opt.verbose
      case 1
        fprintf('%d\tmax: %f\ttarget:%f\tobserved %f\r', iter, max(Yt), target_val, yt(end));
    end
    if opt.plot
        clf
        hold on
        [Z, IZ] = sort(Xs(:,1));
        plot(Z,mu(IZ),'k');
        target_plot = min(mu) + (target-min(target))*(max(mu)-min(mu))/(max(target)-min(target));
        plot(Z,target_plot(IZ),'g');
        plot(Xt(:,1),Yt,'+r');
        drawnow
    end
end

switch opt.verbose
  case 1
    fprintf('\n');
end

end