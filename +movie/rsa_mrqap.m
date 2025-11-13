function stats = rsa_mrqap(S_neural, S_behav, S_cov, varargin)
% RSA_MRQAP  MRQAP for ISC-RSA with optional covariates (FL/DSP), rich output.
% (Memory-safe: uses QR residualization; never forms P = X(X'X)^{-1}X'.)
%
%   Ordinary least squares (dyad-long) + node-label permutations that preserve
%   dyadic dependence (optionally within strata). Supports time-resolved Y
%   (N×N×T), max-T correction across time, rank-based regression, and per-
%   predictor permutation p-values (main-only by default, or all predictors).
%
% USAGE
%   % (1) Mantel-equivalent MRQAP (no covariates)
%   stats = rsa_mrqap(S_neural, S_behav, []);
%
%   % (2) One covariate
%   stats = rsa_mrqap(S_neural, S_behav, S_cov1);
%
%   % (3) Multiple covariates (cell)
%   stats = rsa_mrqap(S_neural, S_behav, {S_cov1, S_cov2}, 'scheme','DSP', 'test','all');
%
% INPUTS
%   S_neural : [N x N] or [N x N x T]  neural similarity (or distance)
%   S_behav  : [N x N]                  behavioral similarity (or distance)
%   S_cov    : [] | {} | [N x N] | cell of [N x N] dyadic covariates
%
% NAME-VALUE PAIRS
%   'scheme'    : 'FL' (default) | 'DSP'         % Freedman–Lane or Double Semi-Partialling
%   'vectorize' : 'lower' (default) | 'upper' | 'alloff'
%   'rank'      : false                          % rank-transform y/x/covs (Spearman-like)
%   'nPerm'     : 5000
%   'twoTailed' : true
%   'maxT'      : true                           % FWER control across time (if T>1)
%   'seed'      : []
%   'strata'    : []                             % length-N stratum labels for restricted perms
%   'expectSym' : true                           % enforce symmetry & zero diag for all matrices
%   'direction' : 'auto'|'similarity'|'distance' % unify sign convention for S_behav vs S_neural
%   'test'      : 'main' (default) | 'all'       % permutation p-values for main only or all preds
%   'names'     : {} or cellstr, length = 1 + 1 + p_cov
%                 e.g., {'(Intercept)','S_behav','SameSex','AgeSim'}
%   'stat'      : 'beta' | 't' (default: 't')    % permutation test statistic
%
% OUTPUT (struct)
%   .beta, .se, .tval               % p×T observed OLS (p = 1+K predictors incl. intercept)
%   .p_main                         % 1×T permutation p for main (S_behav)
%   .p_all                          % p×T permutation p for ALL predictors (if test='all')
%   .permDist_main                  % permutation distribution for main (shape depends on maxT)
%   .permDist_all                   % cell per predictor (only if test='all')
%   .ols                            % struct: n,p,df,sigma,R2,adjR2,SSE,AIC,BIC (per-time arrays)
%   .vecMask                        % logical [N x N] mask of dyads actually analyzed
%   .scheme, .names, .info          % bookkeeping (stat used, rank, nPerm, etc.)
%
% NOTES & GUARANTEES
%   • The SAME dyad set (vecMask) is used in all fits and permutations.
%   • Node-label permutations (rows/cols jointly) preserve dyadic dependence/exchangeability.
%   • Optional stratified permutations (site/batch) & max-T for time-resolved multiplicity.
%   • You may choose β or t as the permutation statistic; t is often more nearly pivotal.


% RSA_MRQAP  MRQAP for ISC-RSA with optional covariates (FL/DSP), rich output.


% -------------------- Parse name-value pairs --------------------
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
ip.addParameter('test','main');                 % 'main' | 'all'
ip.addParameter('names',{});                    % optional column names
ip.addParameter('stat','t', @(s) any(strcmpi(s,{'beta','t','tval'})));
ip.parse(varargin{:});
opt = ip.Results;
if ~isempty(opt.seed), rng(opt.seed); end
useT = any(strcmpi(opt.stat,{'t','tval'}));

% -------------------- Basic checks --------------------
[N1,N2,maybeT] = size(S_neural);
assert(N1==N2, 'S_neural must be square.');
assert(all(size(S_behav,1:2)==[N1 N2]), 'S_behav size mismatch.');
isTime = (ndims(S_neural)==3 && maybeT>1);
T = max(1, maybeT);

% -------------------- Normalize covariate container --------------------
if isempty(S_cov)
    covMats = {};
elseif iscell(S_cov)
    covMats = S_cov;
else
    assert(ismatrix(S_cov) && all(size(S_cov)==[N1 N2]), 'S_cov must be [N x N].');
    covMats = {S_cov};
end
p_cov = numel(covMats);

% -------------------- Names for predictors (incl. intercept) --------------------
if isempty(opt.names)
    names = cell(1, 2 + p_cov);   % [Intercept, S_behav, cov1..covK]
    names{1}='(Intercept)'; names{2}='S_behav';
    for k=1:p_cov, names{2+k}=sprintf('cov%d',k); end
else
    names = opt.names;
    assert(numel(names)==2+p_cov, 'names must be length 1+1+p_cov (incl. intercept).');
end

% -------------------- Symmetry & zero diagonals --------------------
enforce_sym0 = @(A) setdiag0(0.5*(A + A.'));
if opt.expectSym
    S_behav = enforce_sym0(S_behav);
    if isTime
        for t=1:T, S_neural(:,:,t) = enforce_sym0(S_neural(:,:,t)); end
    else
        S_neural = enforce_sym0(S_neural);
    end
    for c=1:p_cov, covMats{c} = enforce_sym0(covMats{c}); end
end

% -------------------- Vectorization (geometry) --------------------
switch lower(string(opt.vectorize))
    case "lower",  geom = tril(true(N1),-1);
    case "upper",  geom = triu(true(N1), 1);
    case "alloff", geom = ~eye(N1);
    otherwise, error('rsa_mrqap:vectorize','Unknown vectorize mode.');
end

% -------------------- Availability mask across ALL variables/T --------------------
avail = isfinite(S_behav);
if isTime
    for t=1:T, avail = avail & isfinite(S_neural(:,:,t)); end
else
    avail = avail & isfinite(S_neural);
end
for c=1:p_cov, avail = avail & isfinite(covMats{c}); end
avail = avail & geom;
if ~any(avail(:)), error('No valid dyads after masking/NA handling.'); end

% -------------------- Vectorize to dyad-long --------------------
if isTime
    y = zeros(nnz(avail), T);
    for t=1:T, tmp = S_neural(:,:,t); y(:,t) = tmp(avail); end
else
    y = S_neural(avail);
end
x_main = S_behav(avail);
Xcov = [];
if p_cov>0
    for c=1:p_cov, Xcov = [Xcov, covMats{c}(avail)]; %#ok<AGROW>
    end
end

% -------------------- Direction alignment (cosmetic) --------------------
if ~strcmpi(opt.direction,'auto')
    if strcmpi(opt.direction,'distance'), x_main = -x_main; end
else
    if isTime, ychk = mean(y,2,'omitnan'); else, ychk = y; end
    rchk = corr(ychk, x_main, 'Type','Spearman', 'Rows','pairwise');
    if ~isnan(rchk) && rchk < 0, x_main = -x_main; end
end

% -------------------- Optional rank transform --------------------
if opt.rank
    x_main = tiedrank(x_main);
    if ~isempty(Xcov)
        for c=1:size(Xcov,2), Xcov(:,c) = tiedrank(Xcov(:,c)); end
    end
    if isTime
        for t=1:T, y(:,t) = tiedrank(y(:,t)); end
    else
        y = tiedrank(y);
    end
end

% -------------------- Designs --------------------
X_full = [ones(numel(x_main),1), x_main, Xcov];   % p_all = 1 + 1 + p_cov
p_all  = size(X_full,2);
useDSP = strcmpi(opt.scheme,'DSP');
if useDSP && isempty(Xcov)
    useDSP = false; % DSP needs ≥1 covariate; fallback to FL
end

% Reduced designs (for FL; per predictor j)
if ~useDSP
    X_red_all = cell(1, p_all);
    for j=1:p_all
        keep = true(1,p_all); keep(j) = false;
        X_red_all{j} = X_full(:, keep);
    end
end

% -------------------- Observed OLS (QR) --------------------
[beta_obs, t_obs, se_obs, ols_info] = ols_full_qr(y, X_full); % p×T
obs_stat = beta_obs; if useT, obs_stat = t_obs; end

stats.beta = beta_obs;
stats.beta_main = beta_obs(2);
stats.se   = se_obs;
stats.tval = t_obs;
stats.ols  = ols_info;

% -------------------- Permutations --------------------
B = opt.nPerm;
pvals_main = NaN(1,T);
pvals_all  = NaN(p_all,T);
permDist_main = [];
permDist_all  = [];

if B>0
    perms_idx = generate_label_perms(N1, B, opt.strata);
    wantAll   = strcmpi(opt.test,'all');

    % ---- MAIN effect p-values (column 2) ----
    if ~useDSP
        % Freedman–Lane for main: reduced model excludes col 2
        X_red_main = X_red_all{2};
        [b_red, ~] = ols_full_qr(y, X_red_main);
        yhat_red = X_red_main * b_red;
        e = y - yhat_red;                                 % (#dyads)×T

        if T>1 && opt.maxT, permDist_main = nan(1,B); else, permDist_main = nan(B,T); end
        for b=1:B
            P   = perms_idx{b};
            e_b = apply_perm_to_residuals(e, N1, avail, P, T);
            y_b = yhat_red + e_b;
            [beta_b, t_b] = ols_full_qr(y_b, X_full);
            stat_b = beta_b(2,:); if useT, stat_b = t_b(2,:); end
            if T>1 && opt.maxT, permDist_main(b) = max_stat(stat_b, opt.twoTailed);
            else,               permDist_main(b,:) = stat_b; end
        end

    else
        % DSP for main: residualize Y and x_main on other predictors (all except col 2) via QR
        keep = true(1,p_all); keep(2)=false;
        Qj = qr_Q(X_full(:,keep));         % thin Q
        y_perp = y  - Qj*(Qj'*y);
        x_perp = X_full(:,2) - Qj*(Qj'*X_full(:,2));

        if T>1 && opt.maxT, permDist_main = nan(1,B); else, permDist_main = nan(B,T); end
        for b=1:B
            P    = perms_idx{b};
            x_b  = apply_perm_to_predictor(x_perp, N1, avail, P);
            [beta_b, t_b] = ols_full_qr(y_perp, [ones(size(x_b)), x_b]);
            stat_b = beta_b(2,:); if useT, stat_b = t_b(2,:); end
            if T>1 && opt.maxT, permDist_main(b) = max_stat(stat_b, opt.twoTailed);
            else,               permDist_main(b,:) = stat_b; end
        end
    end

    % p for main
    if T>1 && opt.maxT
        pvals_main = arrayfun(@(st) perm_pvalue_maxT(permDist_main, st, opt.twoTailed), obs_stat(2,:));
    else
        pvals_main = arrayfun(@(tix) perm_pvalue_vec(permDist_main(:,tix), obs_stat(2,tix), opt.twoTailed), 1:T);
    end

    % ---- ALL predictors p-values (optional) ----
    if wantAll
        permDist_all = cell(1, p_all);
        pvals_all    = NaN(p_all, T);

        for j=1:p_all
            if ~useDSP
                % FL for predictor j
                X_red_j = X_red_all{j};
                [b_red_j, ~] = ols_full_qr(y, X_red_j);
                yhat_red_j = X_red_j * b_red_j;
                e_j = y - yhat_red_j;

                if T>1 && opt.maxT, pd = nan(1,B); else, pd = nan(B,T); end
                for b=1:B
                    P   = perms_idx{b};
                    e_b = apply_perm_to_residuals(e_j, N1, avail, P, T);
                    y_b = yhat_red_j + e_b;
                    [beta_b, t_b] = ols_full_qr(y_b, X_full);
                    stat_b = beta_b(j,:); if useT, stat_b = t_b(j,:); end
                    if T>1 && opt.maxT, pd(b) = max_stat(stat_b, opt.twoTailed);
                    else,               pd(b,:) = stat_b; end
                end
                permDist_all{j} = pd;

                if T>1 && opt.maxT
                    pvals_all(j,:) = arrayfun(@(st) perm_pvalue_maxT(pd, st, opt.twoTailed), obs_stat(j,:));
                else
                    for t=1:T
                        pvals_all(j,t) = perm_pvalue_vec(pd(:,t), obs_stat(j,t), opt.twoTailed);
                    end
                end

            else
                % DSP for predictor j: residualize on all other columns via QR (no P)
                keep = true(1,p_all); keep(j)=false;
                Qj = qr_Q(X_full(:,keep));
                y_perp = y - Qj*(Qj'*y);
                x_perp = X_full(:,j) - Qj*(Qj'*X_full(:,j));

                if T>1 && opt.maxT, pd = nan(1,B); else, pd = nan(B,T); end
                for b=1:B
                    P    = perms_idx{b};
                    x_b  = apply_perm_to_predictor(x_perp, N1, avail, P);
                    [beta_b, t_b] = ols_full_qr(y_perp, [ones(size(x_b)), x_b]);
                    stat_b = beta_b(2,:); if useT, stat_b = t_b(2,:); end
                    if T>1 && opt.maxT, pd(b) = max_stat(stat_b, opt.twoTailed);
                    else,               pd(b,:) = stat_b; end
                end
                permDist_all{j} = pd;

                if T>1 && opt.maxT
                    pvals_all(j,:) = arrayfun(@(st) perm_pvalue_maxT(pd, st, opt.twoTailed), obs_stat(j,:));
                else
                    for t=1:T
                        pvals_all(j,t) = perm_pvalue_vec(pd(:,t), obs_stat(j,t), opt.twoTailed);
                    end
                end
            end
        end
    end
end

% -------------------- Package output --------------------
stats.p_main        = pvals_main;
stats.permDist_main = permDist_main;

if B>0 && strcmpi(opt.test,'all')
    stats.p_all        = pvals_all;
    stats.permDist_all = permDist_all;
else
    stats.p_all        = [];
    stats.permDist_all = [];
end

stats.vecMask  = avail;
stats.scheme   = iff(useDSP,'DSP','FL');
stats.names    = names;
stats.info     = struct('rank',opt.rank,'nPerm',B,'twoTailed',opt.twoTailed,...
                        'maxT',logical(opt.maxT),'strata',opt.strata,'vectorize',opt.vectorize,...
                        'direction',opt.direction,'timeResolved',isTime,...
                        'test',opt.test,'p_all_cols',p_all,'n_dyads',nnz(avail),...
                        'statUsed', iff(useT,'t','beta'));
end % ===== END MAIN =====

% ==================== helpers ====================
function A = setdiag0(A), A(1:size(A,1)+1:end) = 0; end

function Q = qr_Q(X)
% Thin QR and return Q only (n x p). Uses economy-size QR.
    [Q,~] = qr(X,0);
end

function [BETA, TVAL, SE, OLS] = ols_full_qr(Y, X)
% OLS via thin-QR: solves without forming (X'X)^{-1}; small p inverse uses R.
% BETA/TVAL/SE: p×T. Diagnostics per time in OLS.
    if isvector(Y), Y = Y(:); end
    [n,p] = size(X);
    Tt    = size(Y,2);

    % Thin QR
    [Q,R] = qr(X,0);              % X = Q (n×p), R (p×p)
    QtY   = Q' * Y;               % p×T
    % Solve R * BETA = QtY  (upper-triangular back-substitution)
    BETA  = R \ QtY;              % p×T

    % Residuals
    E     = Y - Q*(QtY);
    dof   = max(n - p, 1);
    s2    = sum(E.^2, 1) ./ dof;  % 1×T

    % (X'X)^{-1} = inv(R)*inv(R')  → use inv(R) since p is small
    Rin   = inv(R);
    covB  = Rin * Rin.';          % p×p
    SE    = sqrt(diag(covB)) * sqrt(s2);  % p×T
    TVAL  = BETA ./ SE;

    SSE   = sum(E.^2, 1);
    TSS   = sum( (Y - mean(Y,1)).^2, 1 );
    R2    = 1 - SSE./TSS;
    adjR2 = 1 - (1-R2) .* ((n-1)./max(n-p,1));
    sigma = sqrt(s2);
    k     = p;
    AICt  = n.*log(SSE./n) + 2*k;
    BICt  = n.*log(SSE./n) + k*log(n);

    OLS = struct('n',n,'p',p,'dof',dof,'sigma',sigma,'R2',R2,...
                 'adjR2',adjR2,'SSE',SSE,'AIC',AICt,'BIC',BICt);
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
            P = (1:N)'; %#ok<AGROW>
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

function m = max_stat(vec, twoT), m = twoT*max(abs(vec)) + (~twoT)*max(vec); end

function p = perm_pvalue_maxT(permStats, obs, twoT)
    if twoT, p = (1 + nnz(permStats >= abs(obs))) / (numel(permStats)+1);
    else,     p = (1 + nnz(permStats >= obs))     / (numel(permStats)+1); end
end

function p = perm_pvalue_vec(vecStats, obs, twoT)
    if twoT, p = (1 + nnz(abs(vecStats) >= abs(obs))) / (numel(vecStats)+1);
    else,     p = (1 + nnz(vecStats >= obs))         / (numel(vecStats)+1); end
end

function out = iff(cond,a,b)
    if cond
        out=a; 
    else
        out=b; 
    end
end