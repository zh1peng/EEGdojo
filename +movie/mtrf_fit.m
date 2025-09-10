function out = mtrf_fit(Y, X, cov, fs, varargin)
%MTRF_FIT  Forward mTRF for one continuous run (+ circular-shift surrogates).
%
%   out = mtrf_fit(Y, X, cov, fs, ...)
%
%   PURPOSE
%     Fit a forward TRF (encoding model) predicting neural response Y from
%     time-aligned predictors X (and optional covariates cov). Model performance
%     is evaluated in-run; chance levels are estimated with circular time-shift
%     surrogates. Two permutation modes:
%       (1) shift-only-X (default): X is shifted while cov is kept fixed.
%       (2) shift-both: [X, cov] are concatenated then jointly shifted.
%
%   INPUTS
%     Y    : [T x K] neural response (single subject/run; K comps allowed).
%            Do NOT pass multiple subjects as cell arrays here.
%     X    : [T x P] predictors of interest (can be []).
%     cov  : [T x C] nuisance covariates (can be []; columns appended to X).
%     fs   : sampling rate (Hz).
%
%   NAME-VALUE OPTIONS (defaults in brackets)
%     'dir'            : +1 forward, -1 backward                         [1]
%     'tmin'           : min lag (ms)                                  [-100]
%     'tmax'           : max lag (ms)                                    [500]
%     'lambda'         : ridge parameter                               [2^8]
%     'zeropad'        : use mTRF zero-padding (logical)                [true]
%     'standardize'    : z-score X, cov, Y (within run)                [true]
%
%     Response as envelope
%     'responseMode'   : 'signal' or 'env'                           ['signal']
%     'envLowpassHz'   : LP cutoff for envelope (Hz; []=no LP)            [4]
%     'envLog'         : log-transform envelope                         [true]
%
%     Permutation null (circular shifts)
%     'doPerm'         : compute null with time shifts                 [false]
%     'nPerm'          : # permutations                                   [200]
%     'minShiftSec'    : minimum absolute shift (s)                        [45]
%     'maxShiftSec'    : maximum shift (s; [] -> auto)                    []
%     'rngSeed'        : RNG seed for reproducibility                     []
%     'permShiftBoth'  : if true, jointly shift [X,cov]; else shift X   [false]
%
%   OUTPUT (struct)
%     .fs, .tmin, .tmax, .dir, .lambda, .zeropad
%     .t                : TRF lag axis (ms)
%     .K, .P, .C        : #components in Y, #X columns, #cov columns
%     .model_full       : mTRF model (mTRFtrain output)
%     .perf_true        : table with columns: comp, r, R2, mse
%     .null_R2_full     : [K x nPerm] (if doPerm)
%     .deltaR2, .zscore, .pval  : chance-corrected stats (if doPerm)
%     .bestComp, .deltaR2_best, .deltaR2_meanK
%     .permShiftMode    : 'shift-only-X' or 'shift-both'
%
%   NOTES
%     • This function assumes a single continuous run. For multi-subject/group
%       evaluation with parallelization, implement a wrapper (e.g., group_mtrf_fit).
%     • Choose minShiftSec comfortably larger than total lag span.
%
%   DEPENDENCIES
%     • mTRF-Toolbox: mTRFtrain, mTRFpredict
%     • Signal Processing Toolbox (if envelope low-pass is used)

% ---------- Guard & parse ----------
if iscell(Y) || iscell(X) || iscell(cov)
    error('mtrf_fit expects numeric arrays per single subject/run. Use a group-level wrapper for multiple subjects.');
end

validateattributes(Y,   {'double','single'}, {'2d','nonempty'});
validateattributes(fs,  {'double','single'}, {'scalar','positive'});

if isempty(X),   X = zeros(size(Y,1),0); end
if isempty(cov), cov = zeros(size(Y,1),0); end

[T1,K] = size(Y);
[T2,P] = size(X);
[T3,C] = size(cov);
assert(T1==T2 && T2==T3, 'Y, X, and cov must have the same #rows (time).');
T = T1; %#ok<NASGU>

p = inputParser;
p.addParameter('dir',            1,     @(x)isscalar(x) && any(x==[-1 1]));
p.addParameter('tmin',         -100,    @(x)isscalar(x) && isnumeric(x));
p.addParameter('tmax',          500,    @(x)isscalar(x) && isnumeric(x));
p.addParameter('lambda',        2^8,    @(x)isscalar(x) && isnumeric(x));
p.addParameter('zeropad',       true,   @islogical);
p.addParameter('standardize',   true,   @islogical);

p.addParameter('responseMode', 'signal', @(s)ischar(s)||isstring(s));
p.addParameter('envLowpassHz',  4,       @(x)isempty(x) || (isscalar(x) && x>0));
p.addParameter('envLog',        true,    @islogical);

p.addParameter('doPerm',        false,   @islogical);
p.addParameter('nPerm',         200,     @(x)isscalar(x) && x>=1);
p.addParameter('minShiftSec',    45,     @(x)isscalar(x) && x>=1);
p.addParameter('maxShiftSec',   [],      @(x)isempty(x) || (isscalar(x) && x>0));
p.addParameter('rngSeed',       [],      @(x)isempty(x) || isscalar(x));
p.addParameter('permShiftBoth', false,   @islogical);

p.parse(varargin{:});
opt = p.Results;

% ---------- Response representation ----------
Yproc = Y;
if strcmpi(opt.responseMode,'env')
    Yproc = abs(hilbert(Yproc));
    if opt.envLog, Yproc = log(max(Yproc, eps)); end
    if ~isempty(opt.envLowpassHz)
        [b,a] = butter(4, opt.envLowpassHz/(fs/2), 'low');
        Yproc = filtfilt(b,a,Yproc);
    end
end

% ---------- Standardize ----------
if opt.standardize
    if ~isempty(X),   X   = zscore(X);   end
    if ~isempty(cov), cov = zscore(cov); end
    Yproc = zscore(Yproc);
end

% ---------- Build true-design & fit ----------
Xtrue = [X, cov];           % covariates always included in the true model
Dir   = opt.dir;
tmin  = opt.tmin;
tmax  = opt.tmax;
lam   = opt.lambda;
zp    = opt.zeropad;

model_full = mTRFtrain(Xtrue, Yproc, fs, Dir, tmin, tmax, lam, 'zeropad', zp);
[pred_full, ~] = mTRFpredict(Xtrue, Yproc, model_full, 'zeropad', zp);
[r_true, R2_true, mse_true] = perf_metrics(Yproc, pred_full);

% ---------- Circular-shift null ----------
null_R2_full = [];
deltaR2 = NaN(K,1);
zscoreK = NaN(K,1);
pvalK   = NaN(K,1);

permMode = ternary(opt.permShiftBoth, 'shift-both', 'shift-only-X');

if opt.doPerm
    nPerm = opt.nPerm;

    lag_window_sec = max(0, (tmax - tmin)/1000);
    minShiftSecEff = max(opt.minShiftSec, ceil(lag_window_sec + 5));
    if isempty(opt.maxShiftSec)
        maxShiftSecEff = max(minShiftSecEff+1, floor(T1/fs) - minShiftSecEff);
    else
        maxShiftSecEff = max(opt.maxShiftSec, minShiftSecEff+1);
    end
    minShift = max(1, round(minShiftSecEff * fs));
    maxShift = max(minShift+1, round(maxShiftSecEff * fs));

    if ~isempty(opt.rngSeed), rng(opt.rngSeed); end
    shift_samps = randi([minShift, maxShift], nPerm, 1);

    null_R2_full = nan(K, nPerm);

    if opt.permShiftBoth
        Xcat = [X, cov]; % shift both together
        for pidx = 1:nPerm
            s = shift_samps(pidx);
            Xnull = circshift(Xcat, s, 1);
            mdl_s = mTRFtrain(Xnull, Yproc, fs, Dir, tmin, tmax, lam, 'zeropad', zp);
            pred_s = mTRFpredict(Xnull, Yproc, mdl_s, 'zeropad', zp);
            [~, R2_s] = perf_metrics(Yproc, pred_s);
            null_R2_full(:, pidx) = R2_s(:);
        end
    else
        % shift X only, keep cov fixed
        for pidx = 1:nPerm
            s = shift_samps(pidx);
            Xs    = circshift(X, s, 1);
            Xnull = [Xs, cov];
            mdl_s = mTRFtrain(Xnull, Yproc, fs, Dir, tmin, tmax, lam, 'zeropad', zp);
            pred_s = mTRFpredict(Xnull, Yproc, mdl_s, 'zeropad', zp);
            [~, R2_s] = perf_metrics(Yproc, pred_s);
            null_R2_full(:, pidx) = R2_s(:);
        end
    end

    mu0 = mean(null_R2_full, 2);
    sd0 = std(null_R2_full, 0, 2) + eps;
    deltaR2 = R2_true(:) - mu0;
    zscoreK = deltaR2 ./ sd0;
    pvalK   = (sum(null_R2_full >= R2_true(:), 2) + 1) / (nPerm + 1);
end

% ---------- Package output ----------
out = struct();
out.fs    = fs;
out.tmin  = tmin;
out.tmax  = tmax;
out.dir   = Dir;
out.lambda= lam;
out.zeropad = zp;
out.t     = model_full.t;

out.K = K;
out.P = P;
out.C = C;

perf_true = table((1:K).', r_true(:), R2_true(:), mse_true(:), ...
    'VariableNames', {'comp','r','R2','mse'});
out.perf_true      = perf_true;
out.model_full     = model_full;

out.null_R2_full   = null_R2_full;
out.deltaR2        = deltaR2;
out.zscore         = zscoreK;
out.pval           = pvalK;

[~, bestIdx]       = max(deltaR2);
out.bestComp       = bestIdx;
out.deltaR2_best   = deltaR2(bestIdx);
out.deltaR2_meanK  = mean(deltaR2);

out.permShiftMode  = permMode;
end

% ================= helpers =================

function [rcol, R2col, msecol] = perf_metrics(Y, Yhat)
% Column-wise Pearson r, R^2 (1 - SSE/SST), and MSE
% Y, Yhat: T x K
T = size(Y,1);
Yc    = Y    - mean(Y,1);
Yhatc = Yhat - mean(Yhat,1);

SSE   = sum((Y - Yhat).^2, 1);
SST   = sum(Yc.^2, 1);
msecol= SSE(:) / T;

R2    = 1 - SSE ./ max(SST, eps);
R2col = R2(:);

num   = sum(Yc .* Yhatc, 1);
den   = sqrt(sum(Yc.^2,1) .* sum(Yhatc.^2,1)) + eps;
r     = num ./ den;
rcol  = r(:);
end

function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
