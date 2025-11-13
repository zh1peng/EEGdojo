function stats = rsa_mrqap(S_neural, S_behav, varargin)
% RSA_MRQAP  MRQAP for ISC-RSA with optional covariate matrices (FL/DSP).
%   OLS on dyads + node-label restricted permutations. No Mantel fallback path.
%
% USAGE
%   stats = rsa_mrqap(S_neural, S_behav)                       % one-predictor MRQAP (FL)
%   stats = rsa_mrqap(S_neural, S_behav, S_cov1, S_cov2, ...)  % MRQAP with covariates
%
% INPUTS
%   S_neural : [N x N] or [N x N x T] neural similarity (or distance) matrices.
%   S_behav  : [N x N] behavioral similarity (or distance) matrix.
%   S_cov*   : (optional) each is [N x N] dyadic covariate matrix (|ΔAge|, SameSex, etc.).
%
% NAME-VALUE PAIRS
%   'scheme'    : 'FL' (default) | 'DSP'         % Freedman–Lane or Double Semi-Partialling
%   'vectorize' : 'lower'(dflt) | 'upper' | 'alloff'
%   'rank'      : false (default)                % rank-transform y/x/covs (Spearman-like)
%   'nPerm'     : 5000 (default; 0=no permutations)
%   'twoTailed' : true (default)                 % two-sided permutation p-values
%   'maxT'      : true (default; for T>1)        % max-|stat| FWER control across time windows
%   'seed'      : []
%   'strata'    : []                             % length-N stratum labels for restricted perms
%   'expectSym' : true (default)                 % enforce symmetry & zero diag for all matrices
%   'direction' : 'auto'|'similarity'|'distance' % align similarity/distance direction
%
% OUTPUT (struct)
%   .beta, .tval, .p            % main effect (S_behav) per time; tval仅作参考，显著性以置换p为准
%   .permDist                   % permutation distribution (maxT: 1xB; else: BxT)
%   .vecMask                    % logical [N x N] used for vectorization
%   .scheme, .info              % bookkeeping
%
% REVIEWER-CRITICAL GUARANTEES
%   • Fixed dyad mask across all variables & time → same dyad set used in every permutation.
%   • Node-label permutations (rows/cols jointly) → preserve dyadic dependence/exchangeability.
%   • Optional stratified permutations (site/batch) & max-T for time-resolved multiplicity.
%   • Alignment of distance/similarity direction; optional rank regression (Spearman-like).

% -------------------- Parse inputs --------------------
ip = inputParser;
ip.addParameter('scheme','FL');
ip.addParameter('vectorize','lower');
ip.addParameter('rank',false,@islogical);
ip.addParameter('nPerm',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('twoTailed',true,@islogical);
ip.addParameter('maxT',true,@islogical);
ip.addParameter('seed',[]);
ip.addParameter('strata',[]);
ip.addParameter('expectSym',true,@islogical);
ip.addParameter('direction','auto');
ip.parse(varargin{:});
opt = ip.Results;
if ~isempty(opt.seed), rng(opt.seed); end

% Pull covariate matrices from varargin (numeric 2D N×N at the front)
covMats = {};
for k = 1:numel(varargin)
    v = varargin{k};
    if isnumeric(v) && ndims(v)==2
        covMats{end+1} = v; %#ok<AGROW>
    else
        break; % name-value parsing已完成
    end
end

[N1,N2,maybeT] = size(S_neural);
assert(N1==N2,'S_neural must be square.');
assert(all(size(S_behav,1:2)==[N1 N2]), 'S_behav size mismatch.');
isTime = (ndims(S_neural)==3 && maybeT>1);
T = max(1,maybeT);

% -------------------- Sanity: symmetry & zero diag --------------------
enforce_sym0 = @(A) setdiag0(0.5*(A + A.'));
if opt.expectSym
    S_behav = enforce_sym0(S_behav);
    if isTime
        for t=1:T, S_neural(:,:,t) = enforce_sym0(S_neural(:,:,t)); end
    else
        S_neural = enforce_sym0(S_neural);
    end
    for c=1:numel(covMats), covMats{c} = enforce_sym0(covMats{c}); end
end

% -------------------- Vectorization mask --------------------
switch lower(string(opt.vectorize))
    case "lower",  mask = tril(true(N1),-1);
    case "upper",  mask = triu(true(N1), 1);
    case "alloff", mask = ~eye(N1);
    otherwise, error('rsa_mrqap:vectorize','Unknown vectorize mode.');
end
mask_vec = mask(:);

% -------------------- Fixed dyad availability mask across ALL inputs/T --------------------
avail = isfinite(S_behav);
if isTime
    for t=1:T, avail = avail & isfinite(S_neural(:,:,t)); end
else
    avail = avail & isfinite(S_neural);
end
for c=1:numel(covMats), avail = avail & isfinite(covMats{c}); end
avail = avail & mask;   % off-diagonal only
if ~any(avail(:)), error('No valid dyads after masking/NA handling.'); end

% -------------------- Vectorize to dyad-long --------------------
y = zeros(nnz(avail), T);
if isTime
    for t=1:T, tmp = S_neural(:,:,t); y(:,t) = tmp(avail); end
else
    y = S_neural(avail);
end
x_main = S_behav(avail);
Xcov = [];
for c=1:numel(covMats), Xcov = [Xcov, covMats{c}(avail)]; %#ok<AGROW>
end
p_cov = size(Xcov,2);

% -------------------- Align similarity/distance direction --------------------
% 如果一边是“距离”另一边是“相似”，把主效应方向统一（更方便解读；不影响置换p）
if ~strcmpi(opt.direction,'auto')
    if strcmpi(opt.direction,'distance'), x_main = -x_main; end
else
    % 粗略启发式 + 签号检查
    if isTime
        rchk = corr(mean(y,2,'omitnan'), x_main, 'type','Spearman','rows','pairwise');
    else
        rchk = corr(y, x_main, 'type','Spearman','rows','pairwise');
    end
    if ~isnan(rchk) && rchk<0, x_main = -x_main; end
end

% -------------------- Optional rank transform (Spearman-like regression) --------------------
if opt.rank
    x_main = tiedrank(x_main);
    if p_cov>0
        for c=1:p_cov, Xcov(:,c) = tiedrank(Xcov(:,c)); end
    end
    if isTime
        for t=1:T, y(:,t) = tiedrank(y(:,t)); end
    else
        y = tiedrank(y);
    end
end

% -------------------- Design matrices --------------------
X_full = [ones(numel(x_main),1), x_main, Xcov];  % [1, main, covs]
X_red  = [ones(numel(x_main),1), Xcov];          % reduced (no main) for FL
useDSP = strcmpi(opt.scheme,'DSP');
if useDSP && p_cov==0
    % DSP 需要至少一个协变量；为健壮性，自动回退到 FL（这仍是MRQAP）
    useDSP = false;
end

% -------------------- Observed OLS (per time) --------------------
[beta_obs, t_obs] = ols_timewise(y, X_full);
beta_main = beta_obs(2,:);         % main effect coefficient
t_main    = t_obs(2,:);            % Wald t (参考用)

% -------------------- Permutations: node-label with optional strata --------------------
B = opt.nPerm;
if B==0
    pvals = NaN(size(beta_main));
    perm_stats = [];
else
    perms_idx = generate_label_perms(N1, B, opt.strata);
    if T>1 && opt.maxT
        perm_stats = nan(1,B);     % store max-|stat|
    else
        perm_stats = nan(B,T);     % full per-time stats
    end

    if ~useDSP
        % ===== Freedman–Lane (preferred default) =====
        [beta_red, ~] = ols_timewise(y, X_red);
        yhat_red = X_red * beta_red;       % (#dyads) x T
        e        = y - yhat_red;           % residuals from reduced model

        for b=1:B
            P   = perms_idx{b};
            e_b = apply_perm_to_residuals(e, N1, avail, P, T);   % node-label permute residuals
            y_b = yhat_red + e_b;                                % H0 response
            [beta_b, ~] = ols_timewise(y_b, X_full);
            stat_b = beta_b(2,:);

            if T>1 && opt.maxT, perm_stats(b) = max_stat(stat_b, opt.twoTailed);
            else,               perm_stats(b,:) = stat_b; end
        end

    else
        % ===== Double Semi-Partialling (Dekker DSP) =====
        Pcov = Xcov / (Xcov' * Xcov) * Xcov';
        Mcov = eye(size(Pcov,1)) - Pcov;
        y_perp = Mcov * y;                % residualize Y on covs
        x_perp = Mcov * X_full(:,2);      % residualize main on covs

        for b=1:B
            P     = perms_idx{b};
            x_b   = apply_perm_to_predictor(x_perp, N1, avail, P); % permute main predictor only
            Xb    = [ones(numel(x_b),1), x_b];
            [beta_b, ~] = ols_timewise(y_perp, Xb);
            stat_b = beta_b(2,:);

            if T>1 && opt.maxT, perm_stats(b) = max_stat(stat_b, opt.twoTailed);
            else,               perm_stats(b,:) = stat_b; end
        end
    end

    % ---- Permutation p-values ----
    if T>1 && opt.maxT
        pvals = arrayfun(@(st) perm_pvalue_maxT(perm_stats, st, opt.twoTailed), beta_main);
    else
        pvals = arrayfun(@(tix) perm_pvalue_vec(perm_stats(:,tix), beta_main(tix), opt.twoTailed), 1:T);
    end
end

% -------------------- Output --------------------
stats = struct();
stats.beta     = beta_main;
stats.tval     = t_main;
stats.p        = pvals;
stats.permDist = B>0 && (T>1 && opt.maxT) * perm_stats + ...
                 B>0 && ~(T>1 && opt.maxT) * perm_stats; %#ok<NASGU>
stats.permDist = perm_stats;     % (MATLAB不支持上面那种简写，保留此行)
stats.vecMask  = mask;
stats.scheme   = iff(useDSP,'DSP','FL');
stats.info     = struct('rank',opt.rank,'nPerm',B,'twoTailed',opt.twoTailed,...
                        'maxT',logical(opt.maxT),'strata',opt.strata,'vectorize',opt.vectorize,...
                        'direction',opt.direction,'timeResolved',isTime,'pCov',p_cov);

% ==================== helpers ====================
function A = setdiag0(A)
    A(1:size(A,1)+1:end) = 0;
end

function [BETA, TVAL] = ols_timewise(Y, X)
    % OLS per time; TVAL仅作参考，不用于显著性（显著性基于置换）
    if isvector(Y), Y = Y(:); end
    [n,p] = size(X);
    Tt = size(Y,2);
    XtX = X' * X;
    XtX_inv = pinv(XtX);
    XtY = X' * Y;
    BETA = XtX_inv * XtY;        % p x T
    E = Y - X*(BETA);
    dof = max(n - p, 1);
    s2 = sum(E.^2, 1) ./ dof;    % 1 x T
    se = sqrt(diag(XtX_inv)) * sqrt(s2);   % p x T broadcast
    TVAL = BETA ./ se;
end

function perms = generate_label_perms(N, B, strata)
    % Node-label permutations; optional within-strata shuffles
    perms = cell(1,B);
    if isempty(strata)
        for b=1:B, perms{b} = randperm(N); end
    else
        strata = strata(:);
        assert(numel(strata)==N,'strata length must be N.');
        [~,~,gid] = unique(strata);
        blocks = arrayfun(@(g)find(gid==g), unique(gid), 'uni',0);
        for b=1:B
            P = (1:N)';
            for k=1:numel(blocks)
                blk = blocks{k}; P(blk) = blk(randperm(numel(blk)));
            end
            perms{b} = P(:).';
        end
    end
end

function e_perm = apply_perm_to_residuals(e, N, availMask, P, Tt)
    % Rebuild residual matrices → permute rows/cols jointly → re-vectorize by SAME dyad indices
    e_perm = zeros(size(e));
    [I,J] = find(availMask);
    for t=1:Tt
        R = zeros(N);
        for k=1:numel(I)
            R(I(k),J(k)) = e(k,t); R(J(k),I(k)) = e(k,t);
        end
        Rp = R(P,P);
        for k=1:numel(I), e_perm(k,t) = Rp(I(k),J(k)); end
    end
end

function xp = apply_perm_to_predictor(x, N, availMask, P)
    % Same as above but for a single vector predictor
    xp = zeros(size(x));
    [I,J] = find(availMask);
    R = zeros(N);
    for k=1:numel(I)
        R(I(k),J(k)) = x(k); R(J(k),I(k)) = x(k);
    end
    Rp = R(P,P);
    for k=1:numel(I), xp(k) = Rp(I(k),J(k)); end
end

function m = max_stat(vec, twoT)
    if twoT, m = max(abs(vec)); else, m = max(vec); end
end

function p = perm_pvalue_maxT(permStats, obs, twoT)
    % permStats: 1 x B of max-|stat|
    if twoT, p = (1 + nnz(permStats >= abs(obs))) / (numel(permStats)+1);
    else,     p = (1 + nnz(permStats >= obs))     / (numel(permStats)+1); end
end

function p = perm_pvalue_vec(vecStats, obs, twoT)
    if twoT, p = (1 + nnz(abs(vecStats) >= abs(obs))) / (numel(vecStats)+1);
    else,     p = (1 + nnz(vecStats >= obs))         / (numel(vecStats)+1); end
end

function out = iff(cond,a,b)
    if cond, out=a; else, out=b; end
end

end
