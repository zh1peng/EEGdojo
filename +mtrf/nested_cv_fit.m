function OUT = nested_cv_fit(Y, X, varargin)
% NESTED_CV_FIT  Model-agnostic nested cross-validation + deploy branch (Y,X order).
%
%   OUT = nested_cv_fit(Y, X, 'Name', Value, ...)
%
% Inputs
%   Y : [N x M] targets (supports multiple outputs)
%   X : [N x P] predictors (time/samples x features)
%
% Key Name-Value options
%   'Mode'           : 'assess' | 'deploy' | 'both'    (default: 'assess')
%   'Model'          : 'ridge'|'lasso'|'elasticnet'|'rf'|'gbrt'|struct  (default: 'ridge')
%   'ParamGrid'      : struct hyperparameter grid (e.g., struct('lambda',2.^(0:2:18)))
%   'SelectMetric'   : tuning metric: 'corr'|'r2'|'mse'|function_handle (default: 'corr')
%   'KOuter','KInner': outer/inner contiguous folds (default: 5, 4)
%   'GuardBand'      : #samples trimmed around outer test folds (>= max|lag| for TRF) (default: 0)
%   'Standardize'    : 'train'|'none' (default: 'train')
%   'FoldIndexOuter' : custom outer-fold IDs (length N) (default: [])
%   'UseParallel'    : true/false (default: false; uses parfor if available)
%   'Verbose'        : 0/1 (default: 1)
%
% What it does
%   Branch 'assess' (nested CV):
%     • z-scores X and Y using OUTER-TRAIN stats (once), applies to inner folds & test.
%     • tunes on one metric (e.g., Fisher-z mean corr), but records many on the OUTER-TEST:
%         - Pearson r (per-output + Fisher-z aggregated overall, computed in z-space)
%         - Spearman rho (per-output + mean, computed in z-space)
%         - R², MSE, MAE, RMSE, NRMSE, EV (computed in original Y units)
%     • returns per-fold panels and mean±SE summary (Fisher-z aggregation for r).
%   Branch 'deploy' (fit on all data):
%     • standardizes on ALL data, inner CV selects hypers; for ridge uses 1-SE rule.
%     • returns model + {muX,sdX,muY,sdY} to apply on new data (and invert scaling).
%
% Notes
%   • X and Y are z-scored using OUTER-TRAIN stats; those stats are reused across
%     all INNER folds (no leakage into the OUTER test).
%   • If tuning metric is 'mse'/'r2', model selection is done in ORIGINAL Y units
%     (we unscale predictions with OUTER-TRAIN muY/sdY before scoring).
%     If 'corr', selection uses z-space.
%   • lasso / tree ensembles need Statistics & Machine Learning Toolbox.
%   • Parallelism: set 'UseParallel', true (needs Parallel Computing Toolbox).
%
% Example (assessment with ridge):
%   OUT = nested_cv_fit(Y,X,'Mode','assess','Model','ridge',...
%         'ParamGrid',struct('lambda',2.^(0:2:18)),'KOuter',5,'KInner',4,...
%         'GuardBand',maxLag,'SelectMetric','corr','UseParallel',true);
%
% Example (deploy ridge on all data):
%   OUT = nested_cv_fit(Y,X,'Mode','deploy','Model','ridge',...
%         'ParamGrid',struct('lambda',2.^(0:2:18)),'KInner',5);

% ---------- Parse & checks ----------
ip = inputParser;
ip.addParameter('Mode','assess');
ip.addParameter('Model','ridge');
ip.addParameter('ParamGrid', struct('lambda', 2.^(0:2:18)));
ip.addParameter('SelectMetric','corr');
ip.addParameter('KOuter',5,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
ip.addParameter('KInner',4,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
ip.addParameter('GuardBand',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('Standardize','train',@(s)ischar(s)||isstring(s));
ip.addParameter('FoldIndexOuter',[],@(v) isempty(v) || (isvector(v)&&numel(v)==size(X,1)));
ip.addParameter('UseParallel',false,@islogical);
ip.addParameter('Verbose',1,@(x)isnumeric(x)&&isscalar(x));
ip.parse(varargin{:});
opt = ip.Results;

[N, ~] = size(X);
if size(Y,1) ~= N, error('Y and X must have same #rows.'); end
M = size(Y,2);

% ---------- Model dispatch ----------
[train_fn, predict_fn, default_grid, model_tag] = get_model_dispatch(opt.Model);
grid = merge_grids(default_grid, opt.ParamGrid);
select_metric_name = metric_label(opt.SelectMetric);

OUT = struct();
OUT.settings      = opt;
OUT.model_tag     = model_tag;
OUT.param_grid    = grid;
OUT.select_metric = select_metric_name;

% ---------- Optional parallel pool ----------
if opt.UseParallel
    try
      p = gcp('nocreate'); if isempty(p), parpool; end %#ok<NASGU>
    catch
      if opt.Verbose, warning('Parallel pool not available. Running serial.'); end
      opt.UseParallel = false;
    end
end

%% ===================== Branch 1: ASSESS (nested CV) =====================
if any(strcmpi(opt.Mode, {'assess','both'}))
    if isempty(opt.FoldIndexOuter)
        fold_ids = contiguous_folds(N, opt.KOuter);
    else
        fold_ids = opt.FoldIndexOuter(:);
        opt.KOuter = max(fold_ids);
    end

    OUT.outer_scores  = nan(opt.KOuter,1);      % tuning-metric on OUTER-TEST
    OUT.outer_params  = cell(opt.KOuter,1);
    OUT.predict_outer = nan(size(Y));           % predictions in z-space if Standardize='train'
    OUT.y_outer       = nan(size(Y));           % Y test in z-space if Standardize='train'
    OUT.panels.by_fold= cell(opt.KOuter,1);

    % parallel-friendly temporary holders
    scores_tmp = nan(opt.KOuter,1);
    params_tmp = cell(opt.KOuter,1);
    pred_tmp   = cell(opt.KOuter,1);
    ytmp_z     = cell(opt.KOuter,1);
    panels_tmp = cell(opt.KOuter,1);

    if ~opt.UseParallel
        for k = 1:opt.KOuter
            [scores_tmp(k), params_tmp{k}, pred_tmp{k}, ytmp_z{k}, panels_tmp{k}] = ...
                run_outer_fold(k, fold_ids, Y, X, opt, grid, train_fn, predict_fn, opt.SelectMetric, select_metric_name);
        end
    else
        parfor k = 1:opt.KOuter
            [scores_tmp(k), params_tmp{k}, pred_tmp{k}, ytmp_z{k}, panels_tmp{k}] = ...
                run_outer_fold(k, fold_ids, Y, X, opt, grid, train_fn, predict_fn, opt.SelectMetric, select_metric_name);
        end
    end

    for k = 1:opt.KOuter
        te_full = (fold_ids==k);
        [~, te_mask] = apply_guard(te_full, opt.GuardBand, N);
        OUT.predict_outer(te_mask,:) = pred_tmp{k}; % z-space (if standardized)
        OUT.y_outer(te_mask,:)       = ytmp_z{k};   % z-space (if standardized)
    end
    OUT.outer_scores   = scores_tmp;
    OUT.outer_params   = params_tmp;
    OUT.panels.by_fold = panels_tmp;
    OUT.panels.summary = summarize_panel(OUT.panels.by_fold);
end

%% ===================== Branch 2: DEPLOY (fit on all) ====================
if any(strcmpi(opt.Mode, {'deploy','both'}))
    if strcmpi(opt.Standardize,'train')
        [Xz_all, muX_all, sdX_all] = zscore_train(X);
        [Yz_all, muY_all, sdY_all] = zscore_train(Y);
    else
        Xz_all = X;   muX_all = zeros(1,size(X,2)); sdX_all = ones(1,size(X,2));
        Yz_all = Y;   muY_all = zeros(1,size(Y,2)); sdY_all = ones(1,size(Y,2));
    end

    inner_all = contiguous_folds(N, opt.KInner);
    sel = inner_cv_search_parallel(Yz_all, Xz_all, inner_all, train_fn, predict_fn, grid, opt.SelectMetric, false, opt.UseParallel, muY_all, sdY_all);

    params_star = sel.params;
    % 1-SE rule for ridge: largest lambda within 1 SE of best
    if strcmpi(model_tag,'ridge')
        mv = sel.grid_means; sv = sel.grid_ses; lv = sel.grid_lambda_vec;
        if ~isempty(lv)
            [~, ib] = max(mv);
            thr = mv(ib) - sv(ib);
            cand = find(mv >= thr);
            [~, imax] = max(lv(cand));
            params_star.lambda = lv(cand(imax));
        end
    end

    model_all = train_fn(Xz_all, Yz_all, params_star);

    OUT.deploy = struct();
    OUT.deploy.model  = model_all;
    OUT.deploy.params = params_star;
    OUT.deploy.preproc = struct('muX',muX_all,'sdX',sdX_all,'muY',muY_all,'sdY',sdY_all);
    OUT.deploy.select_metric = select_metric_name;
end
end

%% ========================= Outer-fold runner ============================
function [score_k, params_k, Yhat_z_k, Yte_z_k, panel_k] = run_outer_fold(k, fold_ids, Y, X, opt, grid, train_fn, predict_fn, select_metric, select_metric_name)
% RUN_OUTER_FOLD  One OUTER fold: standardize (train stats), inner-CV tune,
% refit, predict OUTER-TEST, compute tuning metric + full panel.

N = size(X,1);
te_full = (fold_ids==k);
[tr_mask, te_mask] = apply_guard(te_full, opt.GuardBand, N);

Xtr = X(tr_mask,:);  Ytr = Y(tr_mask,:);
Xte = X(te_mask,:);  Yte = Y(te_mask,:);     % Yte in ORIGINAL units (save)

% Standardize once using OUTER-TRAIN stats; reuse across INNER folds & test
if strcmpi(opt.Standardize,'train')
    [Xtr_z, muX, sdX] = zscore_train(Xtr);
    Xte_z = zscore_apply(Xte, muX, sdX);
    [Ytr_z, muY, sdY] = zscore_train(Ytr);
    Yte_z = zscore_apply(Yte, muY, sdY);
else
    Xtr_z = Xtr; Xte_z = Xte; Ytr_z = Ytr; Yte_z = Yte;
    muY = zeros(1,size(Y,2)); sdY = ones(1,size(Y,2)); % for unscale helper
end

% INNER CV on pre-standardized OUTER-TRAIN (no further scaling inside)
inner_ids = contiguous_folds(sum(tr_mask), opt.KInner);
best = inner_cv_search_parallel(Ytr_z, Xtr_z, inner_ids, train_fn, predict_fn, grid, select_metric, false, false, muY, sdY);

% Refit on all OUTER-TRAIN (z-space), test once
model = train_fn(Xtr_z, Ytr_z, best.params);
Yhat_z = predict_fn(model, Xte_z);

% Tuning metric on OUTER-TEST:
%   - 'corr'/'spearman' uses z-space (Yte_z, Yhat_z)
%   - 'mse'/'r2' uses ORIGINAL units (unscale Yhat_z -> Yhat)
score_k = eval_select_metric(select_metric, Yte, Yte_z, Yhat_z, muY, sdY);

% Full panel for OUTER-TEST:
%   - corr/spearman in z-space
%   - losses/R²/EV in original units
Yhat = unscale_from_z(Yhat_z, muY, sdY);
panel_k = metrics_panel_dual(Yte, Yhat, Yte_z, Yhat_z);

params_k = best.params;
Yhat_z_k = Yhat_z;
Yte_z_k  = Yte_z;

if opt.Verbose
    fprintf('[assess] Fold %d | %s=%.4f | params=%s\n', ...
        k, select_metric_name, score_k, jsonencode(params_k));
end
end

%% ====================== Inner CV (parallel aware) =======================
function best = inner_cv_search_parallel(Ytr_z, Xtr_z, inner_ids, train_fn, predict_fn, grid, select_metric, doNotUse, usePar, muY, sdY) %#ok<INUSD>
% INNER_CV_SEARCH_PARALLEL  Grid search across INNER folds on z-standardized data.
% We DO NOT re-standardize inside inner folds; we’re already in OUTER-TRAIN z-space.

plist = grid2list(grid);
K = max(inner_ids);
P = numel(plist);
scores = nan(P, K);

if usePar
    parfor p = 1:P
        scores(p,:) = inner_eval_one_param(Ytr_z, Xtr_z, inner_ids, train_fn, predict_fn, plist{p}, select_metric, muY, sdY);
    end
else
    for p = 1:P
        scores(p,:) = inner_eval_one_param(Ytr_z, Xtr_z, inner_ids, train_fn, predict_fn, plist{p}, select_metric, muY, sdY);
    end
end

mean_scores = mean(scores,2,'omitnan');
[~, idx] = max(mean_scores);

best.params     = plist{idx};
best.cv_scores  = scores(idx,:);
best.mean       = mean_scores(idx);
best.se         = std(scores(idx,:),'omitnan') / sqrt(K);
best.grid_means = mean_scores;
best.grid_ses   = zeros(P,1);
for p=1:P, best.grid_ses(p) = std(scores(p,:),'omitnan') / sqrt(K); end
if isfield(grid,'lambda'), best.grid_lambda_vec = grid.lambda(:); else, best.grid_lambda_vec = []; end
end

function sc = inner_eval_one_param(Ytr_z, Xtr_z, inner_ids, train_fn, predict_fn, params, select_metric, muY, sdY)
% INNER_EVAL_ONE_PARAM  Score one param set across INNER folds (z-space data).
K = max(inner_ids);
sc = nan(1,K);
for k = 1:K
    te_full = (inner_ids==k);
    [tr_mask, te_mask] = apply_guard(te_full, 0, size(Xtr_z,1)); %#ok<ASGLU>
    Xtr2 = Xtr_z(tr_mask,:);  Ytr2 = Ytr_z(tr_mask,:);
    Xte2 = Xtr_z(te_mask,:);  Yte2 = Ytr_z(te_mask,:);

    model = train_fn(Xtr2, Ytr2, params);
    Yhat2_z = predict_fn(model, Xte2);

    % INNER selection metric: corr/spearman in z-space; mse/r2 in ORIGINAL units
    % (Here we don't have the original Yte; but we can unscale with OUTER-TRAIN muY/sdY)
    if ischar(select_metric) || isstring(select_metric)
        sc(k) = eval_select_metric(select_metric, [], Yte2, Yhat2_z, muY, sdY);
    else
        % custom: assume it expects z-space; user can wrap otherwise
        sc(k) = select_metric(Yte2, Yhat2_z);
    end
end
end

%% =========================== Metrics & utils ============================
function score = eval_select_metric(select_metric, Yte_orig, Yte_z, Yhat_z, muY, sdY)
% EVAL_SELECT_METRIC  Apply the right space for selection metric.
switch lower(char(select_metric))
    case {'corr','spearman'}
        % Use z-space (Fisher-z mean corr or mean Spearman)
        if strcmpi(select_metric,'corr')
            score = fisher_mean_colwise_corr(Yte_z, Yhat_z);
        else
            score = spearman_mean(Yte_z, Yhat_z);
        end
    case {'mse','r2'}
        % Use ORIGINAL units: unscale predictions then score
        Yhat = unscale_from_z(Yhat_z, muY, sdY);
        if strcmpi(select_metric,'mse')
            score = -mean((Yte_orig - Yhat).^2, 'all', 'omitnan'); % higher=better
        else
            score = mean(r2_cols(Yte_orig, Yhat), 'omitnan');
        end
    otherwise
        % assume custom handle provided upstream (already applied)
        error('Unknown SelectMetric in eval_select_metric dispatch.');
end
end

function M = metrics_panel_dual(Y_orig, Yhat_orig, Y_z, Yhat_z)
% METRICS_PANEL_DUAL  Record many metrics with correct spaces:
%   - Pearson r / Spearman: z-space (Y_z, Yhat_z)
%   - R², MSE, MAE, RMSE, NRMSE, EV: original units (Y_orig, Yhat_orig)

% Pearson r per output (z-space)
m = size(Y_z,2);
r = nan(1,m);
for j=1:m
    y  = Y_z(:,j) - mean(Y_z(:,j),'omitnan');
    yh = Yhat_z(:,j) - mean(Yhat_z(:,j),'omitnan');
    r(j) = sum(y.*yh,'omitnan') / (sqrt(sum(y.^2,'omitnan')*sum(yh.^2,'omitnan')) + eps);
end
z = atanh(max(min(r,0.999999),-0.999999));
r_overall = tanh(mean(z,'omitnan'));

% Spearman (z-space)
rho = nan(1,m);
for j=1:m
    rho(j) = corr(Y_z(:,j), Yhat_z(:,j), 'Type','Spearman', 'Rows','pairwise');
end
rho_mean = mean(rho,'omitnan');

% Losses & R² (original units)
SSE = sum((Y_orig - Yhat_orig).^2,1);
SST = sum((Y_orig - mean(Y_orig,1,'omitnan')).^2,1) + eps;
R2  = 1 - SSE ./ SST;
R2_mean = mean(R2,'omitnan');

MSE  = mean((Y_orig - Yhat_orig).^2, 'all', 'omitnan');
RMSE = sqrt(MSE);
MAE  = mean(abs(Y_orig - Yhat_orig), 'all', 'omitnan');

sdY = std(Y_orig,0,1,'omitnan'); sdY(sdY==0)=1;
NRMSE = mean( sqrt(mean((Y_orig - Yhat_orig).^2,1,'omitnan')) ./ sdY , 'omitnan');

EV_mean = mean( 1 - var(Y_orig - Yhat_orig,0,1,'omitnan') ./ (var(Y_orig,0,1,'omitnan') + eps), 'omitnan');

M = struct();
M.r_per_output        = r;
M.r_fisher_mean       = mean(z,'omitnan'); % fold-wise z-mean
M.r_overall           = r_overall;
M.spearman_per_output = rho;
M.spearman_mean       = rho_mean;
M.R2_per_output       = R2(:);
M.R2_mean             = R2_mean;
M.MSE                 = MSE;
M.MAE                 = MAE;
M.RMSE                = RMSE;
M.NRMSE               = NRMSE;
M.EV_mean             = EV_mean;
end

function S = summarize_panel(panel_cells)
% SUMMARIZE_PANEL  Aggregate per-fold panels into mean±SE (z-agg for r).
K = numel(panel_cells);
pull = @(f) cellfun(@(s) s.(f), panel_cells, 'uni', 1);

zvec = pull('r_fisher_mean');
S.r_overall_mean = tanh(mean(zvec));
S.r_overall_se   = std(zvec)/sqrt(K);

fields = {'spearman_mean','R2_mean','MSE','MAE','RMSE','NRMSE','EV_mean'};
for i=1:numel(fields)
    v = pull(fields{i});
    S.([fields{i} '_mean']) = mean(v,'omitnan');
    S.([fields{i} '_se'])   = std(v,'omitnan')/sqrt(K);
end
end

%% =========================== Standardization ===========================
function [Xz, mu, sd] = zscore_train(X)
% ZSCORE_TRAIN  Z-score data and return train-time parameters.
%   [Xz, mu, sd] = zscore_train(X)
mu = mean(X, 1, 'omitnan');
sd = std(X, 0, 1, 'omitnan');
sd(sd < eps) = 1;
Xz = (X - mu) ./ sd;
end

function Xz = zscore_apply(X, mu, sd)
% ZSCORE_APPLY  Apply a precomputed z-score transform (train mu/sd).
%   Xz = zscore_apply(X, mu, sd)
if size(sd,2) ~= size(mu,2), error('mu and sd must have same #columns.'); end
sd(sd < eps) = 1;
Xz = (X - mu) ./ sd;
end

function Y = unscale_from_z(Yz, muY, sdY)
% UNSCALE_FROM_Z  Invert z-score to original units using train mu/sd.
Y = Yz .* sdY + muY;
end

%% =========================== Metric helpers ============================
function zmean = fisher_mean_colwise_corr(Yz, Yhat_z)
% FISHER_MEAN_COLWISE_CORR  Fisher-z mean of Pearson r across outputs (z-space).
M = size(Yz,2);
r = zeros(1,M);
for m=1:M
    yc = Yz(:,m)-mean(Yz(:,m),'omitnan');
    hc = Yhat_z(:,m)-mean(Yhat_z(:,m),'omitnan');
    r(m) = sum(yc.*hc,'omitnan') / (sqrt(sum(yc.^2,'omitnan')*sum(hc.^2,'omitnan')) + eps);
end
zmean = mean(atanh(max(min(r,0.999999),-0.999999)),'omitnan');
end

function s = spearman_mean(Yz, Yhat_z)
% SPEARMAN_MEAN  Mean Spearman rho across outputs (z-space).
M = size(Yz,2);
rho = nan(1,M);
for m=1:M
    rho(m) = corr(Yz(:,m), Yhat_z(:,m), 'Type','Spearman', 'Rows','pairwise');
end
s = mean(rho,'omitnan');
end

function r2 = r2_cols(Y, Yhat)
% R2_COLS  Column-wise R² in original units.
SSE = sum((Y-Yhat).^2,1);
SST = sum((Y-mean(Y,1,'omitnan')).^2,1) + eps;
r2  = 1 - SSE./SST;
r2  = r2(:);
end

%% =============================== CV utils ==============================
function ids = contiguous_folds(N, K)
% CONTIGUOUS_FOLDS  Evenly split 1..N into K contiguous folds (no shuffling).
edges = round(linspace(0,N,K+1));
ids = zeros(N,1);
for k=1:K, ids(edges(k)+1:edges(k+1)) = k; end
end

function [tr_mask, te_mask] = apply_guard(te_mask_full, guard, N)
% APPLY_GUARD  Train/test masks with guard band around the test region.
tr_mask = ~te_mask_full;
te_idx = find(te_mask_full);
if isempty(te_idx), te_mask=false(N,1); return; end
i1 = te_idx(1); i2 = te_idx(end);
lo = max(1, i1-guard); hi = min(N, i2+guard);
tr_mask(lo:i1-1) = false; tr_mask(i2+1:hi) = false;
te_mask = false(N,1);
lo2 = min(i1+guard, i2); hi2 = max(i2-guard, i1);
if lo2 <= hi2, te_mask(lo2:hi2) = true; end
end

function L = grid2list(grid)
% GRID2LIST  Cartesian product of grid fields -> list of param structs.
fn = fieldnames(grid);
if isempty(fn), L={struct()}; return; end
vals = cellfun(@(f) grid.(f), fn, 'uni', 0);
[vals{:}] = ndgrid(vals{:}); n = numel(vals{1});
L = cell(n,1);
for i=1:n
    s = struct();
    for j=1:numel(fn)
        v = vals{j}; s.(fn{j}) = v(i);
    end
    L{i} = s;
end
end

function S = merge_grids(A,B)
% MERGE_GRIDS  User grid overrides defaults.
S = A; f = fieldnames(B);
for i=1:numel(f), S.(f{i}) = B.(f{i}); end
end

%% =========================== Model dispatch ============================
function [train_fn, predict_fn, default_grid, tag] = get_model_dispatch(model)
% GET_MODEL_DISPATCH  Return train/predict fns & default grid for a model.
if ischar(model) || isstring(model)
    switch lower(string(model))
        case "ridge"
            tag = 'ridge';
            default_grid = struct('lambda', 2.^(0:2:18));
            train_fn   = @(Xz,Yz,P) train_ridge(Xz, Yz, getOr(P,'lambda',256));
            predict_fn = @(M,Xz) Xz * M.W;

        case "lasso"
            tag = 'lasso'; default_grid = struct('Alpha',1,'Lambda',logspace(-4,2,20));
            train_fn   = @(Xz,Yz,P) train_lasso(Xz,Yz,P);
            predict_fn = @(M,Xz) Xz * M.W;

        case "elasticnet"
            tag = 'elasticnet';
            default_grid = struct('Alpha', linspace(0.1,0.9,5), 'Lambda', logspace(-4,2,18));
            train_fn   = @(Xz,Yz,P) train_lasso(Xz,Yz,P);
            predict_fn = @(M,Xz) Xz * M.W;

        case "rf"
            tag = 'rf';
            default_grid = struct('NumLearningCycles',[100 200], 'MinLeafSize',[1 5 10], 'NumPredictorsToSample',[]);
            train_fn   = @(Xz,Yz,P) train_ens(Xz,Yz,P,'Bag');
            predict_fn = @(M,Xz) predict_ens(M, Xz);

        case "gbrt"
            tag = 'gbrt';
            default_grid = struct('NumLearningCycles',[100 200], 'LearnRate',[0.05 0.1 0.2], 'MinLeafSize',[1 5 10]);
            train_fn   = @(Xz,Yz,P) train_ens(Xz,Yz,P,'LSBoost');
            predict_fn = @(M,Xz) predict_ens(M, Xz);

        otherwise
            error('Unknown Model "%s".', model);
    end
elseif isstruct(model) && isfield(model,'train_fn') && isfield(model,'predict_fn')
    tag = 'custom'; default_grid = struct();
    train_fn = model.train_fn; predict_fn = model.predict_fn;
else
    error('Model must be a known string or a struct with train_fn/predict_fn.');
end
end

%% ======================= Ridge train/predict ============================
function M = train_ridge(Xz, Yz, lambda)
% TRAIN_RIDGE  Fit ridge regression (z-space). Swap in your ridge if desired.
%   M.W maps Xz -> Yz
[n, p] = size(X);
if size(Y,1) ~= n, error('X and Y must have same #rows.'); end
if isvector(Y), Y = reshape(Y, n, []); end
A = X.'*X + lambda*eye(p); A = (A + A.')*0.5;
B = X.'*Y;
try
    R = chol(A); W = R \ (R.' \ B);
catch
    [U,S,V] = svd(X,'econ'); s = diag(S);
    W = V * ( (diag( s ./ (s.^2 + lambda) )) * (U.' * Y) );
end
M = struct('W', W, 'lambda', lambda);
end



%% ====================== Other model trainers/predictors =================
function M = train_lasso(X, Y, P)
% TRAIN_LASSO  Per-output lasso/elastic-net via lasso() (Toolbox).
alpha = getOr(P,'Alpha',1); lam = getOr(P,'Lambda',1e-2);
W = zeros(size(X,2), size(Y,2));
for m=1:size(Y,2)
    [B,~] = lasso(X, Y(:,m), 'Alpha', alpha, 'Lambda', lam);
    W(:,m) = B;
end
M = struct('W', W, 'Params', P);
end

function M = train_ens(X, Y, P, method)
% TRAIN_ENS  Per-output tree ensembles via fitrensemble (Toolbox).
args = {};
if isfield(P,'NumLearningCycles') && ~isempty(P.NumLearningCycles)
    args = [args, {'NumLearningCycles', P.NumLearningCycles}];
end
templ = templateTree('MinLeafSize', getOr(P,'MinLeafSize',1));
if isfield(P,'NumPredictorsToSample') && ~isempty(P.NumPredictorsToSample)
    templ = templateTree('MinLeafSize', getOr(P,'MinLeafSize',1), ...
                         'NumVariablesToSample', P.NumPredictorsToSample);
end
Mcell = cell(1,size(Y,2));
for m=1:size(Y,2)
    Mcell{m} = fitrensemble(X, Y(:,m), 'Method', method, 'Learners', templ, args{:});
end
M = struct('Ensembles', {Mcell}, 'Method', method, 'Params', P);
end

function Yhat = predict_ens(M, X)
% PREDICT_ENS  Predict with per-output ensembles.
Yhat = zeros(size(X,1), numel(M.Ensembles));
for m=1:numel(M.Ensembles), Yhat(:,m) = predict(M.Ensembles{m}, X); end
end

%% ============================== Small utils =============================
function v = getOr(S, f, d)
% GETOR  Get S.(f) or default d if missing/empty.
if isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end
