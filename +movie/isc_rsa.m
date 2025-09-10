function stats = isc_rsa(S_neural, S_behav, varargin)
% ISC_RSA  Intersubject RSA: correlate neural ISC with behavioral similarity.
%
%   stats = isc_rsa(S_neural, S_behav, 'Name',Value,...)
%
% INPUTS
%   S_neural : [N x N] (static) or [N x N x T] (time-resolved) neural similarity/ISC
%   S_behav  : [N x N] behavioral similarity (static)
%
% NAME-VALUE
%   'vectorize' : 'lower' (default) | 'upper' | 'alloff'
%   'method'    : 'spearman' (default) | 'pearson' | 'kendall'
%   'nPerm'     : 5000 (default)
%   'twoTailed' : true (default)
%   'maxT'      : true (default; time-resolved only)
%   'seed'      : [] or numeric
%
% OUTPUT
%   stats.r        : scalar (static) or [1 x T]
%   stats.p        : p-value(s) (max-T corrected if applicable)
%   stats.permDist : permutation distribution (max-|r| if maxT=true)
%   stats.vecMask  : logical mask used for vectorization
%   stats.method, .vectorize, .nPerm, .twoTailed, .maxT

ip = inputParser;
ip.addParameter('vectorize','lower');
ip.addParameter('method','spearman');
ip.addParameter('nPerm',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('twoTailed',true,@islogical);
ip.addParameter('maxT',true,@islogical);
ip.addParameter('seed',[]);
ip.parse(varargin{:});
opt = ip.Results;

if ~isempty(opt.seed), rng(opt.seed); end

[N1,N2,maybeT] = size(S_neural);
assert(N1==N2, 'S_neural must be square');
assert(all(size(S_behav,1:2)==[N1 N2]), 'S_behav size mismatch');

isTime = (ndims(S_neural)==3 && maybeT>1);
T = max(1, maybeT);

switch lower(string(opt.vectorize))
    case "lower",  mask = tril(true(N1),-1);
    case "upper",  mask = triu(true(N1), 1);
    case "alloff", mask = ~eye(N1);
    otherwise, error('isc_rsa:vectorize','Unknown vectorize mode.');
end

corrfun = @(x,y) corr(x(:), y(:), 'type', opt.method, 'rows','complete');

y = S_behav(mask);
r_obs = nan(1,T);

if ~isTime
    x = S_neural(mask);
    r_obs = corrfun(x,y);
else
    for t = 1:T
        x = S_neural(:,:,t);
        r_obs(t) = corrfun(x(mask), y);
    end
end

% Permutation test: shuffle subject labels jointly across rows/cols
nP = opt.nPerm;
if ~isTime
    perm_stats = nan(1,nP);
    for p = 1:nP
        idx = randperm(N1);
        xn = S_neural(idx, idx);
        perm_stats(p) = corrfun(xn(mask), y);
    end
    pval = opt.twoTailed ...
        * (1 + nnz(abs(perm_stats) >= abs(r_obs))) / (nP+1) ...
        + (~opt.twoTailed) * (1 + nnz(perm_stats >= r_obs)) / (nP+1);
    stats = struct('r',r_obs,'p',pval,'permDist',perm_stats, ...
        'vecMask',mask,'method',opt.method,'vectorize',opt.vectorize, ...
        'nPerm',nP,'twoTailed',opt.twoTailed,'maxT',false);

else
    if opt.maxT
        perm_stats = nan(1,nP);
        for p = 1:nP
            idx = randperm(N1);
            rtmp = nan(1,T);
            for t = 1:T
                xn = S_neural(idx, idx, t);
                rtmp(t) = corrfun(xn(mask), y);
            end
            perm_stats(p) = max(abs(rtmp));
        end
        pvals = arrayfun(@(rt) (1 + nnz(perm_stats >= abs(rt))) / (nP+1), r_obs);
    else
        perm_all = nan(nP,T);
        for p = 1:nP
            idx = randperm(N1);
            for t = 1:T
                xn = S_neural(idx, idx, t);
                perm_all(p,t) = corrfun(xn(mask), y);
            end
        end
        pvals = arrayfun(@(t) ...
            (1 + nnz(abs(perm_all(:,t)) >= abs(r_obs(t))))/(nP+1), 1:T);
        perm_stats = perm_all;
    end
    stats = struct('r',r_obs,'p',pvals,'permDist',perm_stats, ...
        'vecMask',mask,'method',opt.method,'vectorize',opt.vectorize, ...
        'nPerm',nP,'twoTailed',opt.twoTailed,'maxT',logical(opt.maxT));
end
end
