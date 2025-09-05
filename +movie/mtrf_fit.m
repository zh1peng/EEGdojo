function out = mtrf_fit(feat_tgt, Y_resp, fs, varargin)
%MTRF_FIT_INDIVIDUAL  Forward mTRF for a single subject + surrogate stats.
%
%   out = mtrf_fit_individual(feat_tgt, Y_resp, fs, ...)
%
%   PURPOSE
%     Fit a forward TRF (encoding model) that predicts neural activity from
%     movie features, using the mTRF-Toolbox. Evaluates performance on the
%     same continuous run, but computes **chance-corrected** effects using
%     **time-shift (circular) surrogates**, which is the standard approach
%     for single-stimulus naturalistic data.
%
%   INPUTS
%     feat_tgt : features (time-synchronized with Y_resp), either:
%                • matrix  [T x P]  (all features; one group called 'all'), or
%                • struct with fields:
%                     .X          [T x P]   feature matrix
%                     .names      {1 x P}   (optional) feature names
%                     .groups     struct    (optional) feature grouping:
%                          .names  {1 x G}  group names
%                          .cols   {1 x G}  cell of index vectors into columns of X
%     Y_resp   : response matrix [T x K]. Typically your SSD outputs:
%                e.g., Y_resp = res12.X_signal(:,1:K) from ssd_pick_at(...).
%     fs       : sampling rate (Hz) for both features and response.
%
%   NAME-VALUE OPTIONS (all optional)
%     'dir'            : model direction, +1 forward (default), -1 backward
%     'tmin'           : minimum lag in ms (default -100)
%     'tmax'           : maximum lag in ms (default  500)
%     'lambda'         : scalar ridge (default 2^8). (Pass a scalar; model is
%                        fit and evaluated on the same run; chance-correction
%                        comes from surrogates rather than CV.)
%     'standardize'    : true/false, z-score X and Y within subject (default true)
%     'responseMode'   : 'signal' (default) or 'env' to use analytic amplitude
%                        of the response (per column) with optional low-pass.
%     'envLowpassHz'   : low-pass cutoff for envelope (default 4), []=no LP
%     'envLog'         : true to log-transform the envelope (default true)
%
%     Surrogates (single-movie significance)
%     'nPerm'          : number of circular-shift surrogates (default 200)
%     'minShiftSec'    : minimum absolute shift (default 45) to break alignment
%     'maxShiftSec'    : max shift (default: duration - 45s)
%     'rngSeed'        : scalar to fix RNG (default [])
%
%     Ablations (unique contributions by feature group)
%     'doAblations'    : true/false (default true)
%
%     Housekeeping
%     'subjectID'      : string or scalar id to store in output (default [])
%
%   OUTPUT (struct) — fields designed for later **group-level** analysis
%     .subjectID         : echo of input
%     .fs, .t, .tmin, .tmax, .lambda, .dir
%     .K, .P, .G         : #components, #features, #groups
%     .groups            : struct with names and column indices used
%     .model_full        : mTRF model struct (weights, metadata) for X_all → Y
%     .perf_true         : table of per-component performance with columns:
%                          comp, r, R2, mse
%     .null_R2           : [K x nPerm]  null R² per component
%     .deltaR2           : [K x 1]      R²_true - mean(null) per component
%     .zscore            : [K x 1]      (R²_true - mean(null)) / std(null)
%     .pval              : [K x 1]      permutation p-value (one-sided)
%
%     If doAblations:
%     .ablations(g)      : struct per group with fields:
%                          name, model, R2_true[Kx1], null_R2[KxnPerm],
%                          deltaR2[Kx1], zscore[Kx1], pval[Kx1]
%
%     Convenience summaries (for quick group stats later):
%     .bestComp          : index of the best component by deltaR2
%     .deltaR2_best      : scalar deltaR2 for best component
%     .deltaR2_meanK     : mean(deltaR2 over K components)
%
%   USAGE EXAMPLE
%     % Y_resp: your top-3 alpha SSD components (T x 3), e.g. res12.X_signal
%     % feat_tgt: T x P matrix or struct as described above
%     out = mtrf_fit_individual(feat_tgt, Y_resp, EEG.srate, ...
%           'tmin',-100,'tmax',500,'lambda',2^8, ...
%           'responseMode','env','envLowpassHz',4,'envLog',true, ...
%           'nPerm',300,'minShiftSec',60,'doAblations',true, ...
%           'subjectID',"S01");
%
%   NOTES
%   • This function **does not** use between-run CV (you have one movie).
%     It evaluates on the same run, and relies on **time-shift surrogates**
%     to estimate chance levels—standard practice for naturalistic data.
%   • Forward mTRF API details are in the official toolbox README and papers. 
%     See mTRFtrain/mTRFpredict/mTRFevaluate.  (Crosse et al., 2016; 2021) 
%     https://github.com/mickcrosse/mTRF-Toolbox
%
%   REFERENCES
%     Crosse MJ et al. (2016) Front Hum Neurosci 10:604 (mTRF toolbox). 
%     Crosse MJ et al. (2021) Front Neurosci 15:705621 (methodological guide).
%     (See toolbox README for links.) 
%
%   © You. BSD-3-Clause dependencies per mTRF-Toolbox.
%
% -------------------------------------------------------------------------

% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('feat_tgt', @(x)isstruct(x) || isnumeric(x));
p.addRequired('Y_resp',   @(x)isnumeric(x) && ismatrix(x));
p.addRequired('fs',       @(x)isnumeric(x) && isscalar(x) && x>0);

p.addParameter('dir',           1, @(x)isscalar(x) && any(x==[-1 1]));
p.addParameter('tmin',       -100, @(x)isscalar(x) && isnumeric(x));
p.addParameter('tmax',        500, @(x)isscalar(x) && isnumeric(x));
p.addParameter('lambda',     2^8,  @(x)isscalar(x) && isnumeric(x));
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

% ---------- Build response (signal vs envelope) ----------
Y = Y_resp;
if strcmpi(opt.responseMode,'env')
    % analytic amplitude per component, optional log and low-pass
    Y = abs(hilbert(Y));
    if opt.envLog, Y = log(max(Y, eps)); end
    if ~isempty(opt.envLowpassHz)
        [b,a] = butter(4, opt.envLowpassHz/(fs/2), 'low');
        Y = filtfilt(b,a,Y);
    end
end

% ---------- Standardize (within subject) ----------
if opt.standardize
    X = zscore(X);
    Y = zscore(Y);
end

% ---------- Set surrogates ----------
nPerm = opt.nPerm;
minShift = round(opt.minShiftSec * fs);
maxShift = isempty(opt.maxShiftSec) * (size(X,1)-minShift) + ...
           ~isempty(opt.maxShiftSec) * round(opt.maxShiftSec * fs);
maxShift = max(minShift+1, maxShift);
if ~isempty(opt.rngSeed), rng(opt.rngSeed); end

% ---------- Fit FULL model on real data ----------
Dir   = opt.dir;
tmin  = opt.tmin;
tmax  = opt.tmax;
lam   = opt.lambda;

model_full = mTRFtrain(X, Y, fs, Dir, tmin, tmax, lam, 'zeropad', 0);
[pred_full, stats_full] = mTRFpredict(X, Y, model_full, 'zeropad', 1);
% mTRFevaluate returns r, R2, mse consistently across K outputs
[R2_true, mse_true] = mTRFevaluate(Y, pred_full);

% ---------- Null distribution via circular-shift (preserve autocorr) ----------
null_R2 = nan(K, nPerm);
for pidx = 1:nPerm
    % choose a random shift in [minShift, maxShift]
    s = randi([minShift, maxShift], 1);
    Xs = circshift(X, s, 1);

    model_s = mTRFtrain(Xs, Y, fs, Dir, tmin, tmax, lam, 'zeropad', 0);
    [pred_s, ~] = mTRFpredict(Xs, Y, model_s, 'zeropad', 1);
    [R2_s, ~] = mTRFevaluate(Y, pred_s);
    null_R2(:, pidx) = R2_s(:);
end

mu0 = mean(null_R2, 2);
sd0 = std(null_R2, 0, 2) + eps;
deltaR2 = R2_true(:) - mu0;
zscoreK = deltaR2 ./ sd0;
pvalK = (sum(null_R2 >= R2_true(:), 2) + 1) / (nPerm + 1); % one-sided

% ---------- Ablations: leave-one-group-out (unique contributions) ----------
ablations = [];
if opt.doAblations
    G = numel(groups.names);
    ablations = repmat(struct('name',[],'model',[],'R2_true',[], ...
        'null_R2',[],'deltaR2',[],'zscore',[],'pval',[]), G, 1);

    for g = 1:G
        cols_keep = setdiff(1:P, groups.cols{g});
        Xg = X(:, cols_keep);

        mg = mTRFtrain(Xg, Y, fs, Dir, tmin, tmax, lam, 'zeropad', 0);
        [pg, ~] = mTRFpredict(Xg, Y, mg, 'zeropad', 1);
        [R2g, ~] = mTRFevaluate(Y, pg);

        null_R2_g = nan(K, nPerm);
        for pidx = 1:nPerm
            s = randi([minShift, maxShift], 1);
            Xsg = circshift(Xg, s, 1);
            mg_s = mTRFtrain(Xsg, Y, fs, Dir, tmin, tmax, lam, 'zeropad', 0);
            [pg_s, ~] = mTRFpredict(Xsg, Y, mg_s, 'zeropad', 1);
            [R2g_s, ~] = mTRFevaluate(Y, pg_s);
            null_R2_g(:, pidx) = R2g_s(:);
        end

        mu0g = mean(null_R2_g, 2);
        sd0g = std(null_R2_g, 0, 2) + eps;
        dR2g = R2g(:) - mu0g;
        zg   = dR2g ./ sd0g;
        pgv  = (sum(null_R2_g >= R2g(:), 2) + 1) / (nPerm + 1);

        ablations(g).name     = groups.names{g};
        ablations(g).model    = mg;
        ablations(g).R2_true  = R2g(:);
        ablations(g).null_R2  = null_R2_g;
        ablations(g).deltaR2  = dR2g;
        ablations(g).zscore   = zg;
        ablations(g).pval     = pgv;
    end
end

% ---------- Package output ----------
out.subjectID = opt.subjectID;
out.fs   = fs;
out.tmin = tmin; out.tmax = tmax; out.dir = Dir; out.lambda = lam;

% time axis in ms from model (mTRF stores lags in model.t)
out.t    = model_full.t; 

out.K = K; out.P = P; out.G = numel(groups.names);
out.groups = groups;

out.model_full = model_full;

perf_true = table((1:K).', stats_full.r(:), R2_true(:), mse_true(:), ...
    'VariableNames', {'comp','r','R2','mse'});
out.perf_true = perf_true;

out.null_R2  = null_R2;
out.deltaR2  = deltaR2;
out.zscore   = zscoreK;
out.pval     = pvalK;

out.ablations = ablations;

[~, bestIdx] = max(deltaR2);
out.bestComp     = bestIdx;
out.deltaR2_best = deltaR2(bestIdx);
out.deltaR2_meanK = mean(deltaR2);

end
