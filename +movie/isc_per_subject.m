function [ISC, meta, ISCpair] = isc_per_subject(Y, varargin)
% ISC_PER_SUBJECT  Inter-subject similarity (2-D only): per-subject LOSO-ISC and pairwise S×S.
%
%   [ISC, meta] = isc_per_subject(Y, ...)
%   [ISC, meta, ISCpair] = isc_per_subject(Y, ...)
%
% INPUT
%   Y : [T x N]  (time x subjects). 2-D only. If you have 3-D (T x D x N),
%       collapse/choose components externally (e.g., average across D) to get T x N.
%
% NAME-VALUE OPTIONS (all optional)
%   'timeResolved' : logical (default: false)
%       If true, compute time-resolved ISC with a sliding window.
%   'winLength'    : positive integer (default: [])
%       Window length (samples) for time-resolved ISC. Required if timeResolved=true.
%       If even is given, it will be incremented by 1 for symmetry.
%   'edgeMode'     : 'shrink' (default) | 'reflect' | 'replicate'
%       Edge handling for time-resolved mode.
%           shrink    = “don’t invent data; accept shorter windows near edges”
%           reflect   = “invent data by mirroring the signal”
%           replicate = “invent data by holding the edge value constant”
%   'pairMetric'   : 'corr' (default) | 'cov'
%       Subject-by-subject similarity metric used **consistently** for LOSO-ISC and ISCpair.
%       - 'corr' → C = corrcoef(windowed Y), bounded [-1,1].
%       - 'cov'  → C = cov(windowed Y), variance-scaled (not bounded).
%   'fisherZ'      : logical (default: true)
%       If true and pairMetric='corr', apply Fisher z-transform **to ISCpair only**.
%       (LOSO-ISC is computed directly from C; do not z-transform C beforehand.)
%   'pairwise'     : logical (default: false)
%       If true and timeResolved=true, compute ISCpair over time (N x N x T).
%       If false in timeResolved mode, ISCpair is [] unless static.
%   'step'         : positive integer (default: 1)
%       Step size (in samples) between successive window centers in time-resolved mode.
%       step = 1 → center at every time point (maximal overlap).
%   'overlap'      : scalar in [0,1) or [] (default: [])
%       Desired overlap fraction between consecutive windows in time-resolved mode.
%       If non-empty, this overrides 'step' via:
%           step = max(1, floor(winLength * (1 - overlap))).
%       Example: winLength=100, overlap=0.8 → step=20 (80% overlap).
%
% OUTPUTS
%   ISC     : If timeResolved=false -> [N x 1]  (per-subject LOSO-like ISC over full time)
%             If timeResolved=true  -> [N x T]  (per-subject ISC at each time center)
%                Only time points used as centers have non-NaN values; others are NaN.
%   meta    : struct with fields:
%               .T, .N, .timeResolved, .winLength, .edgeMode,
%               .pairMetric, .fisherZ, .pairwiseComputed,
%               .step, .overlap, .centerIdx
%   ISCpair : Pairwise subject-by-subject similarity
%             If timeResolved=false -> [N x N]
%             If timeResolved=true  -> [N x N x T] (NaN at unused centers)
%
% NOTES
%   - LOSO-ISC is computed from the N×N similarity matrix C via:
%       Rb_i = 2*sum_{j≠i} C(i,j)
%       Rw_i = (N-1)*C(i,i) + sum_{j≠i} C(j,j)
%       ISC_i = Rb_i / (Rw_i + eps)
%     If pairMetric='corr' (so C has unit diagonal), this simplifies to the
%     **arithmetic mean of off-diagonal correlations** for subject i.
%   - For RSA, use ISCpair (N×N or N×N×T). If averaging correlations, prefer
%     Fisher z on ISCpair slices, average in z, then tanh back.
%
% Fisher-z handling (important):
%   • ISCpair (pairwise):
%       If pairMetric='corr' AND fisherZ=true, ISCpair is returned in
%       Fisher-z space (z = atanh(r)). Averaging and statistical tests on
%       pairwise correlations should be performed in z-space, with results
%       back-transformed using tanh for reporting.
%   • ISC (per-subject LOSO):
%       Always computed directly from C via the ratio Rb/Rw; NO Fisher-z
%       is applied inside this function so ISC remains in its native units
%       (correlation if 'corr', covariance if 'cov').
%       — With pairMetric='corr', ISC_i reduces EXACTLY to the arithmetic
%         mean of off-diagonal correlations for subject i (in r-units).
%       — With pairMetric='cov', ISC is variance-scaled and WILL NOT equal
%         the mean of off-diagonal covariances (by design).
%
% Recommended usage:
%   • For RSA and any averaging over ISCpair: set fisherZ=true, perform
%     all averaging/stats in z-space, then tanh back for display.
%   • For plots and descriptive results: show ISC in correlation units
%     (if 'corr'). If you run parametric stats across subjects on ISC,
%     you MAY z-transform ISC externally (zISC = atanh(rISC)) and report
%     back-transformed estimates.
% v2.1 (2025) – unified metric control; static + time-resolved.
% v2.2 (2025) – avoid to give ISCpair.
% v2.3 (2025) – support time-resolved step / overlap control.

% ---------- Validate input ----------
if ndims(Y) ~= 2
    error('isc_per_subject:2Donly', ...
        'Y must be 2-D [T x N]. If 3-D, collapse to 2-D before calling.');
end
[T, N] = size(Y);

% ---------- Parse options ----------
P = inputParser;
P.addParameter('timeResolved', false, @(x)islogical(x)&&isscalar(x));
P.addParameter('winLength',    [],    @(x)isempty(x) || (isscalar(x) && x>=2 && x==floor(x)));
P.addParameter('edgeMode',     'shrink', @(s)ischar(s) || isstring(s));
P.addParameter('pairMetric',   'corr', @(s) any(strcmpi(string(s),["corr","cov"])));
P.addParameter('fisherZ',      true, @(x)islogical(x)&&isscalar(x));
P.addParameter('pairwise',     false, @(x) isempty(x) || (islogical(x) && isscalar(x)));
P.addParameter('step',         1, @(x)isscalar(x) && x>=1 && x==floor(x));
P.addParameter('overlap',      [], @(x)isempty(x) || (isscalar(x) && x>=0 && x<1));
P.parse(varargin{:});
opt = P.Results;

opt.edgeMode   = lower(string(opt.edgeMode));
opt.pairMetric = lower(string(opt.pairMetric));

% Decide whether to compute pairwise matrices
% Default:
%   - static: compute pairwise if caller asked for 3rd output
%   - time-resolved: do NOT compute pairwise unless explicitly requested
if isempty(opt.pairwise)
    wantPair = (nargout >= 3) && ~opt.timeResolved;
else
    wantPair = (nargout >= 3) && opt.pairwise;
end

meta = struct('T',T,'N',N, ...
              'timeResolved',logical(opt.timeResolved), ...
              'winLength',opt.winLength, ...
              'edgeMode',char(opt.edgeMode), ...
              'pairMetric',char(opt.pairMetric), ...
              'fisherZ',logical(opt.fisherZ), ...
              'pairwiseComputed',logical(wantPair), ...
              'step',[], ...
              'overlap',[], ...
              'centerIdx',[]);

% ---------- Static branch ----------
if ~opt.timeResolved
    % LOSO-ISC from across-subject similarity over full time (consistent metric)
    C = safe_gram(Y, opt.pairMetric);   % N x N (corr or cov)
    ISC = isc_from_C(C);                % [N x 1]

    % Pairwise matrix (same metric) — gated
    if wantPair
        ISCpair = C;
        if opt.fisherZ && opt.pairMetric == "corr"
            ISCpair = atanh(max(min(ISCpair, 0.999999), -0.999999)); % Fisher z
        end
    else
        ISCpair = [];
    end

    % meta.step/overlap irrelevant in static mode
    return;
end

% ---------- Time-resolved branch ----------
if isempty(opt.winLength)
    error('isc_per_subject:winLengthRequired', ...
          'When timeResolved=true, you must provide ''winLength'' (>=2).');
end
w = opt.winLength;
if mod(w,2)==0
    warning('isc_per_subject:evenWindow', ...
        'winLength=%d is even; incremented to %d for symmetric centering.', w, w+1);
    w = w + 1;
end
h = floor((w-1)/2);

% Determine step from overlap or step parameter
if ~isempty(opt.overlap)
    if opt.overlap < 0 || opt.overlap >= 1
        error('isc_per_subject:badOverlap', ...
            'overlap must be in [0,1). Got %.3f.', opt.overlap);
    end
    step = max(1, floor(w * (1 - opt.overlap)));
else
    step = opt.step;
end

% Store in meta (time-resolved)
meta.winLength = w;
meta.step      = step;
meta.overlap   = opt.overlap;
centerIdx      = 1:step:T;
meta.centerIdx = centerIdx(:).';

ISC     = nan(N, T, 'like', Y);      % per-subject ISC
if wantPair
    ISCpair = nan(N, N, T, 'like', Y); % pairwise matrices
else
    ISCpair = [];                      % not computed / not allocated
end

switch opt.edgeMode
    case "shrink"
        % variable window length near edges
        for t = centerIdx
            t1 = max(1, t - h);
            t2 = min(T, t + (w-1 - (t - t1)));
            W = Y(t1:t2, :);                      % [win x N]
            if size(W,1) < 2, continue; end
            C = safe_gram(W, opt.pairMetric);     % N x N
            ISC(:, t) = isc_from_C(C);            % [N x 1]

            if wantPair
                P = C;                             % pairwise
                if opt.fisherZ && opt.pairMetric == "corr"
                    P = atanh(max(min(P, 0.999999), -0.999999));
                end
                ISCpair(:,:,t) = P;
            end
        end

    case {"reflect","replicate"}
        Ypad = pad_time(Y, h, opt.edgeMode);      % [T+2h x N]
        for t = centerIdx
            idx = t:(t+w-1);
            W = Ypad(idx, :);                     % [w x N]
            C = safe_gram(W, opt.pairMetric);
            ISC(:, t) = isc_from_C(C);

            if wantPair
                P = C;
                if opt.fisherZ && opt.pairMetric == "corr"
                    P = atanh(max(min(P, 0.999999), -0.999999));
                end
                ISCpair(:,:,t) = P;
            end
        end

    otherwise
        error('isc_per_subject:edgeMode', 'Unknown edgeMode: %s', opt.edgeMode);
end

end % isc_per_subject


% ==================== Helpers ====================

function C = safe_gram(W, metric)
% Subject-by-subject similarity (corr or cov), robust to short windows.
% W: [T_win x N] -> C: [N x N]
if size(W,1) < 2
    C = nan(size(W,2));
    return;
end
switch string(metric)
    case "corr"
        C = corrcoef(W);   % symmetric, diag=1
    case "cov"
        C = cov(W);        % symmetric, diag=variances
    otherwise
        error('safe_gram:metric','Unknown pairMetric: %s', metric);
end
end

function isc_vec = isc_from_C(C)
% Per-subject LOSO-ISC from an N x N similarity matrix C (corr or cov).
% Rb_i = 2*sum_{j≠i} C(i,j)
% Rw_i = (N-1)*C(i,i) + sum_{j≠i} C(j,j)
% ISC_i = Rb_i / (Rw_i + eps)
N = size(C,1);
diagC = diag(C);
sumOff = sum(C,2) - diagC;            % N x 1
sumDiagOthers = sum(diagC) - diagC;   % N x 1
Rb = 2 .* sumOff;
Rw = (N-1).*diagC + sumDiagOthers;
isc_vec = Rb ./ (Rw + eps);
end

function Ypad = pad_time(Y, h, mode)
% Pad along time by h samples at both ends.
% Y: [T x N] -> Ypad: [T+2h x N]
[T, ~] = size(Y);
if h <= 0, Ypad = Y; return; end
switch string(mode)
    case "reflect"
        hL = min(h, T-1);
        hR = min(h, T-1);
        left  = Y(hL:-1:1, :);
        right = Y(T:-1:(T-hR+1), :);
        Ypad  = [left; Y; right];
    case "replicate"
        left  = repmat(Y(1,:), [h, 1]);
        right = repmat(Y(end,:), [h, 1]);
        Ypad  = [left; Y; right];
    otherwise
        error('pad_time:mode','Unknown pad mode: %s', mode);
end
end
