function out = mtrf_fit(feat_tgt, Y_resp, fs, varargin)
%MTRF_FIT  Forward mTRF for a single continuous run + circular-shift surrogates.
%
%   out = mtrf_fit(feat_tgt, Y_resp, fs, ...)
%
%   PURPOSE
%     Fit a forward temporal response function (TRF; encoding model) that predicts
%     neural activity Y from time-aligned movie features X using the mTRF-Toolbox.
%     Because you typically have one long continuous movie per subject, model
%     performance is evaluated on the same run and **chance levels** are estimated
%     with **circular time-shift surrogates** of X (preserving autocorrelation).
%     Optional, paired **leave-one-group-out ablations** quantify the **unique**
%     contribution of feature groups to prediction performance.
%
%   KEY FEATURES
%     • Forward model: X → Y (lag window [tmin, tmax] ms; ridge regularization).
%     • Optional response envelope: analytic amplitude (+ log/low-pass).
%     • Within-subject standardization (z-score) of X and Y.
%     • Surrogate null via circular shifts of X (paired & reused across tests).
%     • Unique effects via paired nested comparison: FULL vs. REDUCED (–group).
%     • Explicit, toolbox-agnostic metrics per component: r, R², MSE.
%
%   INPUTS
%     feat_tgt : Movie features time-locked to Y (same T rows), either:
%                • numeric [T x P]  — feature matrix (all features = one group).
%                • struct with fields:
%                    .X       [T x P]   feature matrix
%                    .names   {1 x P}   (optional) feature names
%                    .groups  struct    (optional) grouping:
%                              .names {1 x G}  group names
%                              .cols  {1 x G}  column indices into X per group
%
%     Y_resp   : Neural response [T x K], e.g., CorrCA/SSD components for this subject.
%
%     fs       : Sampling rate in Hz (applies to both X and Y).
%
%   NAME-VALUE OPTIONS (all optional; defaults in brackets)
%     'dir'            : Model direction ( +1 = forward, -1 = backward )  [1]
%     'tmin'           : Min lag in ms (negative allowed)                [-100]
%     'tmax'           : Max lag in ms                                     [500]
%     'lambda'         : Ridge parameter (scalar)                          [2^8]
%     'zeropad'        : mTRF zero-padding flag for train + predict (0/1)   [0]
%     'standardize'    : Z-score X and Y within subject (true/false)     [true]
%
%     'responseMode'   : 'signal' (raw Y) or 'env' (analytic envelope)  ['signal']
%     'envLowpassHz'   : Envelope low-pass cutoff Hz ([], no LP)            [4]
%     'envLog'         : Log-transform envelope (true/false)              [true]
%
%     Surrogates (chance level for a single movie)
%     'nPerm'          : # circular-shift surrogates (≥ 20)                [200]
%     'minShiftSec'    : Min absolute shift in seconds                      [45]
%     'maxShiftSec'    : Max shift in seconds ([] → auto)                  [[]]
%                        (Function ensures minShiftSec ≫ lag span and chooses
%                         sensible defaults to avoid trivial alignment.)
%     'rngSeed'        : Scalar RNG seed for reproducibility               [[]]
%
%     Ablations (unique contributions via full vs. reduced comparison)
%     'doAblations'    : Run leave-one-group-out paired ablations         [true]
%
%     Housekeeping
%     'subjectID'      : Subject identifier to echo in output              [[]]
%
%   OUTPUT (struct)
%     .subjectID        : Echo of input
%     .fs               : Sampling rate (Hz)
%     .tmin, .tmax      : Lag window used (ms)
%     .dir, .lambda     : Model direction (+1 forward) and ridge used
%     .zeropad          : Zero-padding flag used for train/predict
%     .t                : Time-lag axis from the fitted model (ms)
%
%     Data/feature summary
%     .K, .P, .G        : #response components, #features, #groups
%     .groups           : Struct with fields .names and .cols actually used
%     .feature_names    : Cellstr feature names if provided (else [])
%
%     Full model (X_all → Y)
%     .model_full       : mTRF model struct returned by mTRFtrain
%     .perf_true        : Table (K rows) with columns:
%                           comp (1..K), r, R2, mse
%     .null_R2_full     : [K x nPerm] null R² per component (circular shifts)
%     .deltaR2          : [K x 1]     R²_true − mean(null) per component
%     .zscore           : [K x 1]     (R²_true − mean(null)) / std(null)
%     .pval             : [K x 1]     One-sided permutation p-value
%
%     Ablations (only if doAblations = true)
%     .ablations(g)     : For each group g:
%                           .name           : group name
%                           .model          : reduced mTRF model (X\group → Y)
%                           .R2_true        : [K x 1]    reduced model R² (true)
%                           .null_R2        : [K x nPerm] reduced null R²
%                           .D_true         : [K x 1]    R²_full − R²_reduced (true)
%                           .D_null         : [K x nPerm]paired null diffs
%                           .deltaR2_unique : [K x 1]    (D_true − mean(D_null))
%                           .zscore         : [K x 1]    (D_true − μ)/σ
%                           .pval           : [K x 1]    One-sided p for unique effect
%
%     Convenience summaries
%     .bestComp         : Index (1..K) of component with max deltaR2
%     .deltaR2_best     : Scalar deltaR2 at best component
%     .deltaR2_meanK    : Mean(deltaR2 over components)
%
%   METHOD OVERVIEW
%     1) (Optional) Transform Y to an amplitude-envelope representation:
%          Y_env = abs(hilbert(Y)); optional log + low-pass.
%     2) Standardize X and Y (z-score across time) within subject.
%     3) Train forward mTRF:  model_full = mTRFtrain(X, Y, fs, dir, tmin, tmax, lambda).
%        Predict on the same run; compute r, R², MSE explicitly per component.
%     4) Null distribution: draw nPerm circular shifts s of X (identical shift set
%        reused across all tests), refit, and recompute R² to obtain null_R2_full.
%     5) Chance-corrected effect sizes: deltaR2, zscore, and one-sided p-values.
%     6) (Optional) Ablation per feature group: train reduced model (X\group),
%        compute **paired** differences D_true = R²_full − R²_reduced and
%        D_null = null_R²_full − null_R²_reduced using the **same shifts**,
%        then report chance-corrected unique contribution.
%
%   USAGE EXAMPLES
%
%     % 1) Minimal usage (all features as one group)
%     out = mtrf_fit(X, Y, fs);
%
%     % 2) Typical setup with lags, ridge, and envelope response
%     out = mtrf_fit(X, Y, fs, ...
%           'tmin', -100, 'tmax', 500, 'lambda', 2^8, ...
%           'responseMode','env', 'envLowpassHz', 4, 'envLog', true, ...
%           'nPerm', 300, 'minShiftSec', 60, 'rngSeed', 42, ...
%           'subjectID', "sub-01");
%
%     % 3) Feature struct with names and groups + ablations
%     F.X = [feat_brightness, feat_motion, feat_loudness, feat_faces]; % T x 4
%     F.names  = {'bright','motion','loud','faces'};
%     F.groups.names = {'visual','auditory','social'};
%     F.groups.cols  = { [1 2], [3], [4] };
%     out = mtrf_fit(F, Y, fs, 'doAblations', true);
%
%     % 4) Backward model (decoder) and consistent zero-padding
%     out = mtrf_fit(X, Y, fs, 'dir', -1, 'zeropad', 0);
%
%     % 5) Reproducible surrogates (fix RNG) and tighter null (more perms)
%     out = mtrf_fit(X, Y, fs, 'rngSeed', 1234, 'nPerm', 1000);
%
%   INTERPRETATION
%     • R²_true is the in-run variance explained by the model per component.
%     • deltaR2 compares R²_true to its circular-shift null mean; zscore and pval
%       quantify evidence that the model predicts above chance.
%     • In ablations, deltaR2_unique reports **how much** performance drops
%       (chance-corrected) when removing a group — a proxy for that group’s
%       unique contribution beyond other features.
%     • If analyzing multiple components (K), consider reporting both bestComp
%       and mean over K; apply multiple-comparison correction if testing per-K p’s.
%
%   BEST PRACTICES
%     • Choose minShiftSec comfortably larger than total lag span (|tmax−tmin|)
%       plus a buffer (this function enforces a safe minimum automatically).
%     • Keep the same preprocessing across subjects/runs for comparability.
%     • Use multivariate X with plausible nuisance features and rely on ablations
%       for **unique** effects. Mass-univariate fits overstate effects when predictors
%       are correlated (common in naturalistic features).
%     • Set rngSeed for reproducibility (esp. if reporting exact p-values).
%
%   DEPENDENCIES
%     • mTRF-Toolbox (Crosse et al., 2016/2021): mTRFtrain, mTRFpredict.
%     • Signal Processing Toolbox for butter/filtfilt (envelope LP) if used.
%
%   REFERENCES
%     Crosse MJ, Di Liberto GM, Bednar A, Lalor EC (2016). The mTRF Toolbox:
%       A MATLAB toolbox for relating neural signals to continuous stimuli.
%       Front. Hum. Neurosci. 10:604.
%     Crosse MJ, Zuk NJ, Di Liberto GM, et al. (2021). Linear modeling of
%       neurophysiological responses to naturalistic stimuli. Front. Neurosci. 15:705621.
%
%   NOTES
%     • Computation scales with nPerm and lag window size; consider parallelizing
%       the permutation loop externally if needed.
%     • This function computes r, R², and MSE explicitly to remain robust
%       across mTRF-Toolbox versions.


% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('feat_tgt', @(x)isstruct(x) || isnumeric(x));
p.addRequired('Y_resp',   @(x)isnumeric(x) && ismatrix(x));
p.addRequired('fs',       @(x)isnumeric(x) && isscalar(x) && x>0);

p.addParameter('dir',           1, @(x)isscalar(x) && any(x==[-1 1]));
p.addParameter('tmin',       -100, @(x)isscalar(x) && isnumeric(x));
p.addParameter('tmax',        500, @(x)isscalar(x) && isnumeric(x));
p.addParameter('lambda',     2^8,  @(x)isscalar(x) && isnumeric(x));
p.addParameter('zeropad',       0, @(x)isscalar(x) && ismember(x,[0 1]));
p.addParameter('standardize', true, @islogical);

p.addParameter('responseMode','signal', @(s)ischar(s)||isstring(s));
p.addParameter('envLowpassHz', 4,  @(x)isempty(x) || (isscalar(x) && x>0));
p.addParameter('envLog',       true, @islogical);

p.addParameter('nPerm',        200, @(x)isscalar(x) && x>=20);
p.addParameter('minShiftSec',   45, @(x)isscalar(x) && x>=1);
p.addParameter('maxShiftSec',  [],  @(x)isempty(x) || (isscalar(x) && x>0));
p.addParameter('rngSeed',      [],  @(x)isempty(x) || isscalar(x));

p.addParameter('doAblations',  true, @islogical);
p.addParameter('subjectID',    [],   @(x)ischar(x)||isstring(x)||isnumeric(x));

p.parse(feat_tgt, Y_resp, fs, varargin{:});
opt = p.Results;

% ---------- Unpack / normalize features ----------
if isstruct(feat_tgt)
    X = feat_tgt.X;
    if isfield(feat_tgt,'names'), feat_names = feat_tgt.names; else, feat_names = []; end
    if isfield(feat_tgt,'groups') && ~isempty(feat_tgt.groups)
        groups = feat_tgt.groups;
        assert(isfield(groups,'cols') && isfield(groups,'names'), ...
            'groups must have .cols (cell of indices) and .names');
    else
        groups.names = {'all'}; groups.cols = {1:size(X,2)};
    end
else
    X = feat_tgt;
    feat_names = [];
    groups.names = {'all'}; groups.cols = {1:size(X,2)};
end
[Tp, P] = size(X);
[Ty, K] = size(Y_resp);
assert(Tp==Ty, 'Feature and response length mismatch.');
T = Tp;

% ---------- Build response (signal vs envelope) ----------
Y = Y_resp;
if strcmpi(opt.responseMode,'env')
    Y = abs(hilbert(Y));                % analytic envelope
    if opt.envLog, Y = log(max(Y, eps)); end
    if ~isempty(opt.envLowpassHz)
        [b,a] = butter(4, opt.envLowpassHz/(fs/2), 'low');
        Y = filtfilt(b,a,Y);
    end
end

% ---------- Standardize (within subject) ----------
if opt.standardize
    X = zscore(X);      % T x P
    Y = zscore(Y);      % T x K
end

% ---------- Set surrogate shifts (paired, reused everywhere) ----------
nPerm = opt.nPerm;
lag_window_sec = max(0, (opt.tmax - opt.tmin)/1000);    % total lag span in seconds
minShiftSecEff = max(opt.minShiftSec, ceil(lag_window_sec + 5));  % ensure well beyond model lags
if isempty(opt.maxShiftSec)
    maxShiftSecEff = max(minShiftSecEff+1, floor(T/fs) - minShiftSecEff);
else
    maxShiftSecEff = max(opt.maxShiftSec, minShiftSecEff+1);
end
minShift = max(1, round(minShiftSecEff * fs));
maxShift = max(minShift+1, round(maxShiftSecEff * fs));
if ~isempty(opt.rngSeed), rng(opt.rngSeed); end
shift_samps = randi([minShift, maxShift], nPerm, 1);    % same shifts reused across all models

% ---------- Fit FULL model on real data ----------
Dir   = opt.dir;
tmin  = opt.tmin;
tmax  = opt.tmax;
lam   = opt.lambda;
zp    = opt.zeropad;

model_full = mTRFtrain(X, Y, fs, Dir, tmin, tmax, lam, 'zeropad', zp);
[pred_full, ~] = mTRFpredict(X, Y, model_full, 'zeropad', zp);

% Explicit metrics (toolbox-agnostic)
[r_true, R2_true, mse_true] = perf_metrics(Y, pred_full);

% ---------- Null distribution via circular-shift (preserve autocorr) ----------
null_R2_full = nan(K, nPerm);
for pidx = 1:nPerm
    s = shift_samps(pidx);
    Xs = circshift(X, s, 1);
    model_s = mTRFtrain(Xs, Y, fs, Dir, tmin, tmax, lam, 'zeropad', zp);
    pred_s  = mTRFpredict(Xs, Y, model_s, 'zeropad', zp);
    [~, R2_s] = perf_metrics(Y, pred_s);
    null_R2_full(:, pidx) = R2_s(:);
end
mu0_full = mean(null_R2_full, 2);
sd0_full = std(null_R2_full, 0, 2) + eps;
deltaR2  = R2_true(:) - mu0_full;
zscoreK  = deltaR2 ./ sd0_full;
pvalK    = (sum(null_R2_full >= R2_true(:), 2) + 1) / (nPerm + 1); % one-sided

% ---------- Ablations: paired (unique contribution) ----------
ablations = [];
if opt.doAblations
    G = numel(groups.names);
    ablations = repmat(struct( ...
        'name',[],'model',[],'R2_true',[], ...
        'null_R2',[],'D_true',[],'D_null',[], ...
        'deltaR2_unique',[],'zscore',[],'pval',[]), G, 1);

    for g = 1:G
        cols_keep = setdiff(1:P, groups.cols{g});
        Xg = X(:, cols_keep);

        % True reduced model
        mg = mTRFtrain(Xg, Y, fs, Dir, tmin, tmax, lam, 'zeropad', zp);
        pred_g = mTRFpredict(Xg, Y, mg, 'zeropad', zp);
        [~, R2g_true] = perf_metrics(Y, pred_g);

        % Paired null: use the very same shifts as for the full model
        null_R2_g = nan(K, nPerm);
        for pidx = 1:nPerm
            s = shift_samps(pidx);
            Xsg = circshift(Xg, s, 1);
            mg_s = mTRFtrain(Xsg, Y, fs, Dir, tmin, tmax, lam, 'zeropad', zp);
            pred_gs = mTRFpredict(Xsg, Y, mg_s, 'zeropad', zp);
            [~, R2g_s] = perf_metrics(Y, pred_gs);
            null_R2_g(:, pidx) = R2g_s(:);
        end

        % Paired unique effect: difference full − reduced
        D_true = R2_true(:) - R2g_true(:);                       % K x 1
        D_null = null_R2_full - null_R2_g;                       % K x nPerm
        muD    = mean(D_null, 2);
        sdD    = std(D_null, 0, 2) + eps;

        % One-sided p: improvement from including group > 0
        pD     = (sum(D_null >= D_true, 2) + 1) / (nPerm + 1);
        zD     = (D_true - muD) ./ sdD;

        ablations(g).name            = groups.names{g};
        ablations(g).model           = mg;
        ablations(g).R2_true         = R2g_true(:);
        ablations(g).null_R2         = null_R2_g;
        ablations(g).D_true          = D_true;
        ablations(g).D_null          = D_null;
        ablations(g).deltaR2_unique  = D_true - muD;    % chance-corrected unique effect
        ablations(g).zscore          = zD;
        ablations(g).pval            = pD;
    end
end

% ---------- Package output ----------
out.subjectID = opt.subjectID;
out.fs   = fs;
out.tmin = tmin; out.tmax = tmax; out.dir = Dir; out.lambda = lam; out.zeropad = zp;

% time axis in ms from model (mTRF stores lags in model.t)
out.t    = model_full.t;

out.K = K; out.P = P; out.G = numel(groups.names);
out.groups = groups;
out.feature_names = feat_names;

out.model_full = model_full;

perf_true = table((1:K).', r_true(:), R2_true(:), mse_true(:), ...
    'VariableNames', {'comp','r','R2','mse'});
out.perf_true = perf_true;

out.null_R2_full = null_R2_full;
out.deltaR2      = deltaR2;
out.zscore       = zscoreK;
out.pval         = pvalK;

out.ablations = ablations;

[~, bestIdx] = max(deltaR2);
out.bestComp       = bestIdx;
out.deltaR2_best   = deltaR2(bestIdx);
out.deltaR2_meanK  = mean(deltaR2);

end % ---- main ----


% ================= helpers =================

function [rcol, R2col, msecol] = perf_metrics(Y, Yhat)
% Column-wise Pearson r, R^2 (1 - SSE/SST), and MSE
% Y, Yhat: T x K
T = size(Y,1);
Yc    = Y   - mean(Y,1);
Yhatc = Yhat- mean(Yhat,1);

SSE = sum((Y - Yhat).^2, 1);
SST = sum(Yc.^2, 1);
msecol = SSE(:) / T;

% Guard: if SST ~ 0 (flat), set R2 = 0
R2 = 1 - SSE ./ max(SST, eps);
R2col = R2(:);

% Pearson r per column
num = sum(Yc .* Yhatc, 1);
den = sqrt(sum(Yc.^2,1) .* sum(Yhatc.^2,1)) + eps;
r = num ./ den;
rcol = r(:);
end
