function [ISC, meta, ISCpair] = isc_per_subject(Y, varargin)
% ISC_PER_SUBJECT  Inter-subject correlation (2-D only): per-subject scores and pairwise S×S.
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
%   'pairMetric'   : 'corr' (default) | 'cov'
%       Pairwise subject-by-subject similarity metric (corr recommended).
%   'fisherZ'      : logical (default: false)
%       If true and pairMetric='corr', apply Fisher z-transform to pairwise matrices.
%
% OUTPUTS
%   ISC     : If timeResolved=false -> [N x 1]  (per-subject ISC over full time)
%             If timeResolved=true  -> [N x T]  (per-subject ISC at each time center)
%   meta    : struct with fields: .T, .N, .timeResolved, .winLength, .edgeMode,
%                                 .pairMetric, .fisherZ
%   ISCpair : Pairwise subject-by-subject similarity
%             If timeResolved=false -> [N x N]
%             If timeResolved=true  -> [N x N x T]
%
% NOTES
%   - Per-subject ISC follows R_b / R_w from the across-subject covariance:
%       C = cov(windowed Y)  -> N x N
%       Rb_i = 2*sum_{j≠i} C(i,j)
%       Rw_i = (N-1)*C(i,i) + sum_{j≠i} C(j,j)
%       ISC_i = Rb_i / (Rw_i + eps)
%   - Pairwise matrices use corrcoef (or cov) within each window.
%   - For time-resolved mode with 'shrink', windows near edges shrink; for
%     'reflect'/'replicate', fixed-length windows are used via padding.
%
% EXAMPLES
%   % Static ISC + pairwise:
%   [ISC, meta, C] = isc_per_subject(Y);           % ISC: [N x 1], C: [N x N]
%
%   % Time-resolved (e.g., 5 s at Fs):
%   [ISCtr, meta, Ctr] = isc_per_subject(Y, 'timeResolved',true, ...
%       'winLength', round(5*Fs), 'edgeMode','reflect');  % ISCtr: [N x T], Ctr: [N x N x T]
%
%   % Fisher-z pairwise (for RSA/regression):
%   [~, ~, Cz] = isc_per_subject(Y, 'fisherZ',true);
%
% v2.0 (2025) – 2D-only design; static + time-resolved per-subject ISC and pairwise matrices.

% ---------- Validate input ----------
if ndims(Y) ~= 2
    error('isc_per_subject:2D only Y must be 2-D [T x N]. If 3-D, collapse to 2-D before calling.');
end
[T, N] = size(Y);

% ---------- Parse options ----------
P = inputParser;
P.addParameter('timeResolved', false, @(x)islogical(x)&&isscalar(x));
P.addParameter('winLength',    [],    @(x)isempty(x) || (isscalar(x) && x>=2 && x==floor(x)));
P.addParameter('edgeMode',     'shrink', @(s)ischar(s) || isstring(s));
P.addParameter('pairMetric',   'corr', @(s) any(strcmpi(string(s),["corr","cov"])));
P.addParameter('fisherZ',      false, @(x)islogical(x)&&isscalar(x));
P.parse(varargin{:});
opt = P.Results;

opt.edgeMode   = lower(string(opt.edgeMode));
opt.pairMetric = lower(string(opt.pairMetric));

meta = struct('T',T,'N',N, ...
              'timeResolved',logical(opt.timeResolved), ...
              'winLength',opt.winLength, ...
              'edgeMode',char(opt.edgeMode), ...
              'pairMetric',char(opt.pairMetric), ...
              'fisherZ',logical(opt.fisherZ));

% ---------- Static branch ----------
if ~opt.timeResolved
    % Per-subject ISC from across-subject covariance over full time
    C = safe_cov(Y);                 % N x N
    ISC = isc_from_cov(C);           % [N x 1]

    % Pairwise matrix
    ISCpair = pairwise_from_timeseries(Y, opt.pairMetric); % [N x N]
    if opt.fisherZ && opt.pairMetric == "corr"
        ISCpair = atanh(max(min(ISCpair, 0.999999), -0.999999)); % Fisher z
    end
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

ISC     = nan(N, T, 'like', Y);   % per-subject ISC
ISCpair = nan(N, N, T, 'like', Y); % pairwise matrices

switch opt.edgeMode
    case "shrink"
        % variable window length near edges
        for t = 1:T
            t1 = max(1, t - h);
            t2 = min(T, t + (w-1 - (t - t1)));
            W = Y(t1:t2, :);               % [win x N]
            if size(W,1) < 2, continue; end
            C = safe_cov(W);               % N x N
            ISC(:, t)   = isc_from_cov(C); % [N x 1]
            P = pairwise_from_timeseries(W, opt.pairMetric);
            if opt.fisherZ && opt.pairMetric == "corr"
                P = atanh(max(min(P, 0.999999), -0.999999));
            end
            ISCpair(:,:,t) = P;
        end

    case {"reflect","replicate"}
        Ypad = pad_time(Y, h, opt.edgeMode); % [T+2h x N]
        for t = 1:T
            idx = t:(t+w-1);
            W = Ypad(idx, :);               % [w x N]
            C = safe_cov(W);
            ISC(:, t)   = isc_from_cov(C);
            P = pairwise_from_timeseries(W, opt.pairMetric);
            if opt.fisherZ && opt.pairMetric == "corr"
                P = atanh(max(min(P, 0.999999), -0.999999));
            end
            ISCpair(:,:,t) = P;
        end

    otherwise
        error('isc_per_subject:edgeMode', 'Unknown edgeMode: %s', opt.edgeMode);
end

end % isc_per_subject

% ==================== Helpers ====================

function C = safe_cov(W)
% Column-wise covariance (subjects), robust to short windows.
% W: [T_win x N] -> C: [N x N]
if size(W,1) < 2
    C = nan(size(W,2));
else
    C = cov(W);
end
end

function isc_vec = isc_from_cov(C)
% Subject-resolved ISC from an N x N covariance matrix.
% Rb_i = 2*sum_{j≠i} C(i,j)
% Rw_i = (N-1)*C(i,i) + sum_{j≠i} C(j,j)
% ISC_i = Rb_i / (Rw_i + eps)
N = size(C,1);
diagC = diag(C);
sumOff = sum(C,2) - diagC;             % N x 1
sumDiagOthers = sum(diagC) - diagC;     % N x 1
Rb = 2 .* sumOff;
Rw = (N-1).*diagC + sumDiagOthers;
isc_vec = Rb ./ (Rw + eps);
end

function P = pairwise_from_timeseries(W, metric)
% Pairwise subject-by-subject similarity within window
% W: [T_win x N] -> P: [N x N]
switch string(metric)
    case "corr"
        if size(W,1) < 2
            P = nan(size(W,2));
        else
            P = corrcoef(W); % symmetric, diag=1
        end
    case "cov"
        P = cov(W);
    otherwise
        error('pairwise_from_timeseries:metric Unknown pairMetric: %s', metric);
end
end

function Ypad = pad_time(Y, h, mode)
% Pad along time by h samples at both ends.
% Y: [T x N] -> Ypad: [T+2h x N]
[T, ~] = size(Y);
if h <= 0, Ypad = Y; return; end
switch mode
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
