function out = ridge_nested_cv(X, y, varargin)
% RIDGE_NESTED_CV  Nested cross-validation for ridge regression with configurable inner metric.
%
%   out = ridge_nested_cv(X, y, 'Kouter', 5, 'Kinner', 4, ...
%                         'purge', 0, 'lambdas', logspace(-3,3,21), ...
%                         'inner_metric', 'corr')
%
% Inputs
%   X : [T x p] predictor matrix
%   y : [T x 1] response vector
%
% Name-Value Pairs
%   'Kouter'       : number of outer folds (default 5)
%   'Kinner'       : number of inner folds (default 4)
%   'purge'        : samples to purge around each outer test block (default 0)
%   'lambdas'      : vector of λ values to test (default logspace(-3,3,21))
%   'inner_metric' : selection metric for inner CV. One of:
%                    'corr' (default; Pearson r on z-space targets),
%                    'mse'  (negative MSE in original units),
%                    'mae'  (negative MAE in original units),
%                    'r2'   (R² in original units),
%                    or a function handle @(y_true,y_pred,y_true_z,y_pred_z)->scalar
%                    where *higher is better*.
%
% Notes
%   • X and y are z-scored using OUTER-TRAIN stats; those stats are used
%     across all inner folds (no leakage into outer test).
%   • If metric is 'mse'/'mae'/'r2', evaluation is in original y units
%     (predictions are unscaled with OUTER-TRAIN stats). For 'corr', evaluation
%     is in z-space (common practice).
%
% Output (struct)
%   .r_folds, .r_mean, .r2_folds, .r2_mean, .mse_folds, .mse_mean
%   .lambda_folds : [Kouter x 1] selected λ per outer fold
%   .yhat_oos     : [T x 1] out-of-sample predictions (original y units)
%
% See also: rcb.blocked_folds, rcb.ridge_regress, rcb.zscore_train, rcb.zscore_apply

ip = inputParser;
ip.addParameter('Kouter', 5, @(x) isnumeric(x) && x>=2);
ip.addParameter('Kinner', 4, @(x) isnumeric(x) && x>=2);
ip.addParameter('purge',  0, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('lambdas', logspace(-3, 3, 21), @isnumeric);
ip.addParameter('inner_metric', 'corr');  % 'corr' | 'mse' | 'mae' | 'r2' | function_handle
ip.parse(varargin{:});
opt = ip.Results;

[T, ~] = size(X);
folds_out = rcb.blocked_folds(T, opt.Kouter, opt.purge);

r_out        = nan(opt.Kouter, 1);
r2_out       = nan(opt.Kouter, 1);
mse_out      = nan(opt.Kouter, 1);
lambda_star  = nan(opt.Kouter, 1);
yhat_all     = nan(T, 1, 'like', y);

% Normalize inner_metric into a scorer function that returns "higher is better"
scorer = make_scorer(opt.inner_metric);

for ko = 1:opt.Kouter
    % ----- Outer split -----
    tr_idx = folds_out(ko).trainIdx;
    te_idx = folds_out(ko).testIdx;

    Xtr = X(tr_idx, :);
    ytr = y(tr_idx, :);
    Xte = X(te_idx, :);
    yte = y(te_idx, :);

    % Z-score X using OUTER-TRAIN stats; apply to test
    [Xtrz, muX, sdX] = rcb.zscore_train(Xtr);
    Xtez = rcb.zscore_apply(Xte, muX, sdX);

    % Z-score y using OUTER-TRAIN stats; apply to test
    [ytrz, muy, sdy] = rcb.zscore_train(ytr);
    ytez = rcb.zscore_apply(yte, muy, sdy);

    % ----- Inner CV for λ selection -----
    folds_in = rcb.blocked_folds(size(Xtrz, 1), opt.Kinner, 0);
    mean_score = nan(numel(opt.lambdas), 1);

    for li = 1:numel(opt.lambdas)
        lam  = opt.lambdas(li);
        scr  = nan(opt.Kinner, 1);

        for ki = 1:opt.Kinner
            tri = folds_in(ki).trainIdx;   % inner-train (subset of outer-train)
            vei = folds_in(ki).testIdx;    % inner-val  (subset of outer-train)

            % Fit on inner-train (z-scored with OUTER-TRAIN stats)
            Wi = rcb.ridge_regress(Xtrz(tri, :), ytrz(tri, :), lam);

            % Predict on inner-val (z-space)
            yv_z = Xtrz(vei, :) * Wi;

            % Prepare targets in both spaces
            yv_true_z = ytrz(vei, :);
            % Unscale predictions/targets to original units
            yv_true   = ytr(vei, :);
            yv        = yv_z * sdy + muy;

            % Score (higher is better)
            scr(ki) = scorer(yv_true, yv, yv_true_z, yv_z);
        end

        mean_score(li) = mean(scr, 'omitnan');
    end

    % Select λ with the best mean inner-CV score (higher is better)
    [~, bestIdx] = max(mean_score);
    lam = opt.lambdas(bestIdx);
    lambda_star(ko) = lam;

    % ----- Refit on full OUTER-TRAIN with best λ; evaluate on OUTER-TEST -----
    W = rcb.ridge_regress(Xtrz, ytrz, lam);

    % Predictions (z-scored target space -> original units)
    yhat_z       = Xtez * W;
    yhat_unscaled= yhat_z * sdy + muy;

    % Store OOS predictions
    yhat_all(te_idx) = yhat_unscaled;

    % Metrics
    r_out(ko)  = corr(yhat_z, ytez, 'rows', 'complete');  % z-space
    ss_res     = sum((yte - yhat_unscaled).^2, 'omitnan');
    ss_tot     = sum((yte - mean(yte, 'omitnan')).^2, 'omitnan');
    r2_out(ko) = 1 - (ss_res / ss_tot);
    mse_out(ko)= mean((yte - yhat_unscaled).^2, 'omitnan');
end

out = struct( ...
    'r_folds',      r_out, ...
    'r_mean',       mean(r_out,  'omitnan'), ...
    'r2_folds',     r2_out, ...
    'r2_mean',      mean(r2_out, 'omitnan'), ...
    'mse_folds',    mse_out, ...
    'mse_mean',     mean(mse_out,'omitnan'), ...
    'lambda_folds', lambda_star, ...
    'yhat_oos',     yhat_all);

end

% ---------- helpers ----------
function scorer = make_scorer(metric)
% Returns a function handle: @(y_true,y_pred,y_true_z,y_pred_z) -> scalar (higher is better)
    if isa(metric, 'function_handle')
        scorer = metric;  % user guarantees "higher is better"
        return;
    end
    validatestring(metric, {'corr','mse','mae','r2'});
    switch lower(metric)
        case 'corr'
            % Pearson r in z-space (common for TRF/encoding)
            scorer = @(y_true,y_pred,y_true_z,y_pred_z) ...
                corr(y_pred_z, y_true_z, 'rows', 'complete');
        case 'mse'
            % Negative MSE in original units (higher is better)
            scorer = @(y_true,y_pred,~,~) ...
                -mean((y_true - y_pred).^2, 'omitnan');
        case 'mae'
            % Negative MAE in original units (higher is better)
            scorer = @(y_true,y_pred,~,~) ...
                -mean(abs(y_true - y_pred), 'omitnan');
        case 'r2'
            % R² in original units (higher is better)
            scorer = @(y_true,y_pred,~,~) ...
                local_r2(y_true, y_pred);
    end
end

function r2 = local_r2(y_true, y_pred)
    ss_res = sum((y_true - y_pred).^2, 'omitnan');
    mu     = mean(y_true, 'omitnan');
    ss_tot = sum((y_true - mu).^2, 'omitnan');
    r2     = 1 - (ss_res / ss_tot);
end
