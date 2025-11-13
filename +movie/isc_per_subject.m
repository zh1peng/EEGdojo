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
%   'pairMetric'   : 'corr' (default) | 'cov'
%       Subject-by-subject similarity metric used **consistently** for LOSO-ISC and ISCPair.
%       - 'corr' → C = corrcoef(windowed Y), bounded [-1,1].
%       - 'cov'  → C = cov(windowed Y), variance-scaled (not bounded).
%   'fisherZ'      : logical (default: true)
%       If true and pairMetric='corr', apply Fisher z-transform **to ISCpair only**.
%       (LOSO-ISC is computed directly from C; do not z-transform C beforehand.)
%
% OUTPUTS
%   ISC     : If timeResolved=false -> [N x 1]  (per-subject LOSO-like ISC over full time)
%             If timeResolved=true  -> [N x T]  (per-subject ISC at each time center)
%   meta    : struct with fields: .T, .N, .timeResolved, .winLength, .edgeMode,
%                                 .pairMetric, .fisherZ
%   ISCpair : Pairwise subject-by-subject similarity
%             If timeResolved=false -> [N x N]
%             If timeResolved=true  -> [N x N x T]
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
% Rationale:
%   Pairwise correlations are not additive; Fisher-z (atanh) makes them
%   approximately normal/additive for averaging. LOSO ISC is commonly
%   reported in correlation units for interpretability and comparability.
%
% Recommended usage:
%   • For RSA and any averaging over ISCPair: set fisherZ=true, perform
%     all averaging/stats in z-space, then tanh back for display.
%   • For plots and descriptive results: show ISC in correlation units
%     (if 'corr'). If you run parametric stats across subjects on ISC,
%     you MAY z-transform ISC externally (zISC = atanh(rISC)) and report
%     back-transformed estimates.
% v2.1 (2025) – unified metric control; static + time-resolved.
% v2.2 (2025) – avoid to give ISCpair.
% ---------- Validate input ----------
if ndims(Y) ~= 2
    error('isc_per_subject:2Donly', 'Y must be 2-D [T x N]. If 3-D, collapse to 2-D before calling.');
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
              'pairwiseComputed',logical(wantPair));

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

ISC     = nan(N, T, 'like', Y);      % per-subject ISC
if wantPair
    ISCpair = nan(N, N, T, 'like', Y); % pairwise matrices
else
    ISCpair = [];                      % not computed / not allocated
end

switch opt.edgeMode
    case "shrink"
        % variable window length near edges
        for t = 1:T
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
                ISCPair(:,:,t) = P;
            end
        end

    case {"reflect","replicate"}
        Ypad = pad_time(Y, h, opt.edgeMode);      % [T+2h x N]
        for t = 1:T
            idx = t:(t+w-1);
            W = Ypad(idx, :);                     % [w x N]
            C = safe_gram(W, opt.pairMetric);
            ISC(:, t) = isc_from_C(C);

            if wantPair
                P = C;
                if opt.fisherZ && opt.pairMetric == "corr"
                    P = atanh(max(min(P, 0.999999), -0.999999));
                end
                ISCPair(:,:,t) = P;
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







%% Sanity checks
% %% ---------- Inputs ----------
% % Y is assumed to be [T x D x N] from CorrCA projections, where you want D=1.
% % If you already have [T x N], set Y2D = Y; and skip the squeeze line.
% Y2D = squeeze(Y(:,1,:));   % [T x N], using component #1

% %% ---------- 1) ISC from CorrCA package (cov-based) ----------
% % CorrCA package function 'isc_per_subject' typically expects [T x N] (or its own format).
% % Replace this call with the correct CorrCA API you use if different.
% ISC1 = isc_per_subject(Y);   % returns [N x 1] LOSO-like ISC (cov-based)
% ISCcov1 = ISC1(:,1);           % ensure it's [N x 1]

% %% ---------- 2) Your unified function: cov and corr ----------
% [ISCcov,  meta_cov,  ISCPair_cov]  = movie.isc_per_subject(Y2D, 'pairMetric','cov');
% [ISCcorr, meta_corr, ISCPair_corr] = movie.isc_per_subject(Y2D, 'pairMetric','corr');

% % Basic shape checks
% assert(isvector(ISCcov1) && isvector(ISCcov) && isvector(ISCcorr), 'ISC outputs must be vectors.');
% assert(numel(ISCcov1)==numel(ISCcov) && numel(ISCcov)==numel(ISCcorr), 'ISC vectors must have same length.');

% N = numel(ISCcov1);

% %% ---------- 3) Numeric sanity: CorrCA(cov) vs ours(cov) ----------
% % They need not be IDENTICAL (different edge/mean conventions), but should be close if formulas match.
% % Report MAE and max|diff|.
% diff_cov = ISCcov - ISCcov1;
% mae_cov  = mean(abs(diff_cov),'omitnan');
% max_cov  = max(abs(diff_cov),[],'omitnan');
% fprintf('CorrCA(cov) vs Ours(cov): MAE=%.4g, max|Δ|=%.4g\n', mae_cov, max_cov);

% %% ---------- 4) Recover LOSO from ISCpair (corr): should match ----------
% % off-diagonal mask (not strictly needed with vectorized sum below)
% N = size(ISCPair_corr,1);

% % Row-wise mean in z-space (exclude diagonal)
% sum_offdiag_z = sum(ISCPair_corr, 2) - diag(ISCPair_corr);   % N x 1
% row_means_z   = sum_offdiag_z ./ (N - 1);                    % N x 1

% % Convert back to correlation units (r) for LOSO comparison
% ISCcorr_from_pair = tanh(row_means_z);                        % N x 1 (r-space)

% % Compare in r-space (LOSO from corr vs row-mean(z)->tanh)
% diff_corr = ISCcorr - ISCcorr_from_pair;
% mae_corr  = mean(abs(diff_corr), 'omitnan');
% max_corr  = max(abs(diff_corr), [], 'omitnan');
% fprintf('[r] LOSO(corr) vs row-mean z(ISCPair)->tanh: MAE=%.4g, max|Δ|=%.4g\n', mae_corr, max_corr);

% % (Optional) also compare in z-space for completeness
% % LOSO(corr) lives in r; convert to z with clipping
% ISCcorr_z = atanh(max(min(ISCcorr, 0.999999), -0.999999));
% diff_corr_z = ISCcorr_z - row_means_z;
% mae_corr_z  = mean(abs(diff_corr_z), 'omitnan');
% max_corr_z  = max(abs(diff_corr_z), [], 'omitnan');
% fprintf('[z] atanh(LOSO(corr)) vs row-mean z(ISCPair): MAE=%.4g, max|Δ|=%.4g\n', mae_corr_z, max_corr_z);
% %% ---------- 5) (Optional) Covariance row-mean (for reference only; not equal to LOSO-cov) ----------
% row_means_cov = (sum(ISCPair_cov - diag(diag(ISCPair_cov)), 2) / (N-1));  % plain mean of off-diag covariances
% diff_cov_mean = ISCcov - row_means_cov;
% mae_cov_mean  = mean(abs(diff_cov_mean),'omitnan');
% max_cov_mean  = max(abs(diff_cov_mean),[],'omitnan');
% fprintf('LOSO(cov) vs row-mean cov (not expected equal): MAE=%.4g, max|Δ|=%.4g\n', mae_cov_mean, max_cov_mean);

% %% ---------- 6) Plots ----------
% % A) Trace comparison (treating subject index as x for static LOSO values)
% figure('Name','ISC (static) comparisons'); 
% tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

% % A1: cov traces (CorrCA vs Ours)
% nexttile; hold on; box on; grid on;
% plot(ISCcov1, 'LineWidth',1.5, 'DisplayName','CorrCA cov');
% plot(ISCcov,  'LineWidth',1.5, 'DisplayName','Ours cov');
% xlabel('Subject'); ylabel('ISC (cov)'); title('Covariance LOSO');
% legend('Location','best'); set(gca,'FontSize',11);

% % A2: corr traces (Ours)
% nexttile; hold on; box on; grid on;
% plot(ISCcorr, 'LineWidth',1.5, 'DisplayName','Ours corr');
% xlabel('Subject'); ylabel('ISC (corr)'); title('Correlation LOSO');
% legend('Location','best'); set(gca,'FontSize',11);

% % A3: differences cov (Ours - CorrCA)
% nexttile; hold on; box on; grid on;
% stem(diff_cov, 'filled', 'DisplayName','Ours - CorrCA');
% xlabel('Subject'); ylabel('\Delta ISC (cov)'); title('Covariance Δ');
% yline(0,'k-'); legend('Location','best'); set(gca,'FontSize',11);

% % A4: LOSO(corr) vs row-mean(Fisher-z ISCPair) scatter
% nexttile; hold on; box on; grid on;
% scatter(ISCcorr_from_pair, ISCcorr, 25, 'filled');
% refline(1,0); % 45° line
% xlabel('Row-mean z(ISCPair) → tanh'); ylabel('LOSO (corr)');
% title(sprintf('Corr equivalence: MAE=%.3g, max|Δ|=%.3g', mae_corr, max_corr));
% set(gca,'FontSize',11);

% %% ---------- 7) Tiny utility: pretty difference plot for corr ----------
% figure('Name','Correlation: LOSO vs from ISCPair'); hold on; box on; grid on;
% plot(diff_corr, 'LineWidth',1.2);
% yline(0,'k-'); 
% xlabel('Subject'); ylabel('\Delta (LOSO - row-mean z(ISCPair)→tanh)');
% title('Corr: LOSO minus Fisher-z row-mean of ISCPair');
% set(gca,'FontSize',12);

