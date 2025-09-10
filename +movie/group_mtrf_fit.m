function TRFres = group_mtrf_fit(Y, X, cov, xNames, covNames, fs, varargin)
% GROUP_MTRF_FIT  Multivariate mTRF per subject using ALL X and ALL cov as predictors.
%
%   TRFres = group_mtrf_fit(Y, X, cov, xNames, covNames, fs, ...)
%
% INPUTS
%   Y        : [T x S] neural response (single run per subject; S subjects)
%   X        : [T x P] predictors of interest (features; can be [])
%   cov      : [T x C] covariates (can be [])
%   xNames   : 1xP cellstr/string; names for X columns (or [])
%   covNames : 1xC cellstr/string; names for cov columns (or [])
%   fs       : sampling rate (Hz)
%
% NAME–VALUE (local & forwarded to mtrf_fit)
%   'useParallel'    : run per-subject fits with parfor            [true]
%
%   Any additional Name–Value pairs are passed to mtrf_fit, e.g.:
%     'dir','tmin','tmax','lambda','zeropad','standardize',...
%     'responseMode','envLowpassHz','envLog',...
%     'doPerm','nPerm','minShiftSec','maxShiftSec','rngSeed','permShiftBoth'
%
% OUTPUT (single struct)
%   .t           : lag axis (ms)
%   .predNames   : names for predictors (X followed by cov)
%   .xIdx        : indices of X predictors in predNames (1..P)
%   .covIdx      : indices of cov predictors (P+1..P+C)
%   .W           : [S x nLags x (P+C)] subject TRFs for all predictors
%   .mu          : [ (P+C) x nLags ] group mean across subjects
%   .se          : [ (P+C) x nLags ] group SE across subjects
%   .std         : [ (P+C) x nLags ] group SD across subjects
%   .perf        : table per subject with: subj, r, R2, mse
%   .permStats   : (if doPerm) table per subject with: subj, deltaR2, z, p
%   .meta        : struct with fs, names, options used
%
% NOTE
%   • Fits one multivariate model per subject with ALL X and ALL cov.
%     To evaluate a specific X column, call this function multiple times
%     with your own X slices (e.g., X(:,i) with full cov).

% ---------------- Parse & validate ----------------
ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('useParallel', true, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

% Convert unmatched options to Name–Value for mtrf_fit
passNV = namedargs2nv(ip.Unmatched);

% Shapes
validateattributes(Y,   {'double','single'}, {'2d','nonempty'});
validateattributes(X,   {'double','single'}, {'2d'});
validateattributes(cov, {'double','single'}, {'2d'});

[T1,S] = size(Y);
[T2,P] = size(X);
[T3,C] = size(cov);
if isempty(X),   P = 0; end
if isempty(cov), C = 0; end

if S==1
    error('group_mtrf_fit:SingleSubject', ...
        'Y is a single time series. Use mtrf_fit(Y, X, cov, fs, ...) for one subject.');
end

% Time-length agreement
if P>0 && T2~=T1
    error('Time length mismatch: Y is %d rows, X is %d rows.', T1, T2);
end
if C>0 && T3~=T1
    error('Time length mismatch: Y is %d rows, cov is %d rows.', T1, T3);
end

% Names
if nargin < 4 || isempty(xNames)
    xNames = compose("X%02d", 1:P);
else
    xNames = cellstr(xNames);
    assert(numel(xNames)==P, 'xNames length (%d) must equal size(X,2)=%d.', numel(xNames), P);
end
if nargin < 5 || isempty(covNames)
    covNames = compose("cov%02d", 1:C);
else
    covNames = cellstr(covNames);
    assert(numel(covNames)==C, 'covNames length (%d) must equal size(cov,2)=%d.', numel(covNames), C);
end

validateattributes(fs, {'double','single'}, {'scalar','positive'});

predNames = [xNames(:).' covNames(:).'];
xIdx   = 1:P;
covIdx = (P+1):(P+C);

% ---------------- Pilot fit to get lag axis & dims ----------------
pilot = movie.mtrf_fit(Y(:,1), X, cov, fs, passNV{:});
t_axis = pilot.t(:).';
nLags  = numel(t_axis);

w0 = pilot.model_full.w;                  % [nPred x nLags]
if size(w0,2) ~= nLags
    error('Unexpected weight shape from mtrf_fit: got [%d x %d], expected nLags=%d in 2nd dim.', size(w0,1), size(w0,2), nLags);
end
nPred = size(w0,1);
if nPred ~= (P + C)
    warning('Predictor count in model (%d) differs from size(X,2)+size(cov,2)=%d. Check inputs.', nPred, P+C);
end

% ---------------- Per-subject fits (parallel) ----------------
W = zeros(S,  nPred, nLags);
r_s = nan(S,1); R2_s = nan(S,1); mse_s = nan(S,1);
dR2_s = nan(S,1); z_s = nan(S,1); p_s = nan(S,1);

usePar = opt.useParallel && license('test','Distrib_Computing_Toolbox');

if usePar
    parfor s = 1:S
        out = movie.mtrf_fit(Y(:,s), X, cov, fs, passNV{:});
        if numel(out.t) ~= nLags || any(out.t(:).' ~= t_axis)
            error('Inconsistent lag axis across subjects; align mTRF options.');
        end
        W(s,:,:) = out.model_full.w;

        r_s(s)   = out.perf_true.r(1);
        R2_s(s)  = out.perf_true.R2(1);
        mse_s(s) = out.perf_true.mse(1);

        if ~isempty(out.deltaR2) && ~all(isnan(out.deltaR2))
            dR2_s(s) = out.deltaR2(1);
            z_s(s)   = out.zscore(1);
            p_s(s)   = out.pval(1);
        end
    end
else
    for s = 1:S
        out = movie.mtrf_fit(Y(:,s), X, cov, fs, passNV{:});
        if numel(out.t) ~= nLags || any(out.t(:).' ~= t_axis)
            error('Inconsistent lag axis across subjects; align mTRF options.');
        end
        W(s,:,:) = out.model_full.w;

        r_s(s)   = out.perf_true.r(1);
        R2_s(s)  = out.perf_true.R2(1);
        mse_s(s) = out.perf_true.mse(1);

        if ~isempty(out.deltaR2) && ~all(isnan(out.deltaR2))
            dR2_s(s) = out.deltaR2(1);
            z_s(s)   = out.zscore(1);
            p_s(s)   = out.pval(1);
        end
    end
end

% ---------------- Group summaries ----------------
mu  = squeeze(mean(W, 1, 'omitnan')).';    % [nPred x nLags]
sd  = squeeze(std(W,  0, 1, 'omitnan')).'; % [nPred x nLags]
nEff= sum(~any(isnan(squeeze(W(:,:,1))),2));
se  = sd ./ max(1, sqrt(nEff));

% ---------------- Package output ----------------
TRFres = struct();
TRFres.t           = t_axis;
TRFres.predNames   = predNames;
TRFres.xIdx        = xIdx;
TRFres.covIdx      = covIdx;
TRFres.W           = W;          % [S x nLags x nPred]
TRFres.mu          = mu;         % [nPred x nLags]
TRFres.se          = se;         % [nPred x nLags]
TRFres.std         = sd;         % [nPred x nLags]

TRFres.perf = table((1:S).', r_s, R2_s, mse_s, ...
    'VariableNames', {'subj','r','R2','mse'});

if any(~isnan(dR2_s))
    TRFres.permStats = table((1:S).', dR2_s, z_s, p_s, ...
        'VariableNames', {'subj','deltaR2','z','p'});
else
    TRFres.permStats = table();
end

TRFres.meta = struct( ...
    'fs', fs, ...
    'xNames', {xNames}, ...
    'covNames', {covNames}, ...
    'options', ip.Unmatched );

end

% -------- helper: struct of NV -> flat NV cell --------
function nv = namedargs2nv(S)
if isempty(S)
    nv = {};
    return;
end
fn = fieldnames(S);
nv = cell(1, numel(fn)*2);
for k = 1:numel(fn)
    nv{2*k-1} = fn{k};
    nv{2*k}   = S.(fn{k});
end
end
