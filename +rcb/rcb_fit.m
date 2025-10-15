function out = rcb_fit(y, Xraw, fs, lags_ms, centers_ms, width_ms, varargin)
% RCB_FIT  Fit a TRF-style encoding model with flexible design choices.
%
%   out = rcb_fit(y, Xraw, fs, lags_ms, centers_ms, width_ms, ...
%                 'DesignType','basis'|'lagged'|'raw', ...
%                 'Kouter',5,'Kinner',4,'PurgeSec',1.0, ...
%                 'Lambdas',logspace(-3,3,21), 'InnerMetric','corr', ...
%                 'AddIntercept',false, 'VariableNames',{})
%
% Notes:
%   • 'basis'  -> raised-cosine basis over lags (reconstruct TRF via B).
%   • 'lagged' -> classic lagged design (TRF is direct reshape of weights).
%   • 'raw'    -> no temporal expansion; W_trf = [].

ip = inputParser;
ip.addParameter('DesignType','basis',@(s) any(strcmpi(s,{'basis','lagged','raw'})));
ip.addParameter('Kouter',      5, @(x) isnumeric(x) && x>=2);
ip.addParameter('Kinner',      4, @(x) isnumeric(x) && x>=2);
ip.addParameter('PurgeSec',  1.0, @(x) isnumeric(x) && isscalar(x) && x>=0);
ip.addParameter('Lambdas', logspace(-3,3,21), @isnumeric);
ip.addParameter('InnerMetric','corr');
ip.addParameter('AddIntercept', false, @(x) islogical(x) || isnumeric(x));
ip.addParameter('VariableNames', {}, @(x) isstring(x) || iscellstr(x));
ip.addParameter('fs', 250, @(x) isnumeric(x) && isscalar(x) && x>=0);
ip.parse(varargin{:});
opt = ip.Results;

% --- 0) Coerce inputs and basic checks ---
if istable(Xraw)
    Xmat = table2array(Xraw);
    featNames = string(Xraw.Properties.VariableNames(:)).';
elseif ismatrix(Xraw)
    Xmat = Xraw;
    F = size(Xmat,2);
    if ~isempty(opt.VariableNames)
        featNames = string(opt.VariableNames(:)).';
        assert(numel(featNames) == F, ...
            'VariableNames (%d) must match Xraw columns (%d).', numel(featNames), F);
    else
        featNames = compose("feat%02d", 1:F);
    end
else
    error('Xraw must be a numeric matrix or a table.');
end
[T1, F] = size(Xmat);
T2 = size(y,1);
assert(T1==T2, 'Xraw and y must have the same number of rows.');

% --- 1) Build design according to DesignType ---
designType = lower(opt.DesignType);
B = [];                 % basis matrix (if used)
L = numel(lags_ms);     % requested lag support length (may be unused for 'raw')

switch designType
    case 'basis'
        % Build raised-cosine basis and convolve features with it
        B = rcb.build_cosine_basis(lags_ms, centers_ms, width_ms);   % [L x K]
        [Lchk, K] = size(B); %#ok<ASGLU>
        assert(Lchk==L, 'length(lags_ms) must equal size(B,1).');
        Xb = rcb.convolve_basis_design(Xmat, B);                     % [T x (F*K)]
    case 'lagged'
        % Classic lagged design: concatenate X shifted by each lag
        % Expect rcb.lagged_design(X, lags_ms) -> [T x (F*L)], with NaNs on edges
        Xb = rcb.lagged_design(Xmat, lags_ms, 'fs', fs, 'units','ms', 'pad','nan');                      % [T x (F*L)]
        K = L;                   % treat as K=L with identity basis conceptually
    case 'raw'
        % Use raw predictors directly (no temporal expansion)
        Xb = Xmat;               % [T x F]
        K = 1;                   % one "temporal slice"
        L = 0;                   % no lag axis
    otherwise
        error('Unknown DesignType: %s', opt.DesignType);
end

% Optional intercept
if opt.AddIntercept
    Xb = [Xb, ones(size(Xb,1),1, 'like', Xb)];
    featNames = [featNames, "intercept"];
end

% --- 2) Handle edge effects / NaNs BEFORE CV ---
bad = any(isnan(Xb), 2) | isnan(y);
if any(bad)
    Xb = Xb(~bad, :);
    y  = y(~bad, :);
    dropped = nnz(bad);
else
    dropped = 0;
end

T = size(Xb,1);
if T < (opt.Kouter * 2)
    warning('After dropping %d invalid rows, only T=%d remain. CV may be unstable.', dropped, T);
end

% --- 3) Nested CV for unbiased performance ---
purge_pts = round(opt.PurgeSec * fs);   % rcb.blocked_folds expects samples
cvres = rcb.ridge_nested_cv(Xb, y, ...
    'Kouter', opt.Kouter, 'Kinner', opt.Kinner, ...
    'purge', purge_pts, 'lambdas', opt.Lambdas, ...
    'inner_metric', opt.InnerMetric);

% --- 4) Final refit on all (cleaned) data with λ = median(λ*)
lam_star = median(cvres.lambda_folds, 'omitnan');
[Xz, muX, sdX] = rcb.zscore_train(Xb);
[yz, muY, sdY] = rcb.zscore_train(y);
W_basis = rcb.ridge_regress(Xz, yz, lam_star);       % column vector or [.. x 1]
yhat_final = (Xz * W_basis) * sdY + muY;

% --- 5) Reconstruct TRF in lag space when meaningful ---
if opt.AddIntercept
    W_core = W_basis(1:end-1, :);   % drop intercept
else
    W_core = W_basis;
end

switch designType
    case 'basis'
        % [(F*K) x 1] -> [K x F] -> [L x F]
        W_trf = B * reshape(W_core, [], F);
    case 'lagged'
        % direct reshape: [(F*L) x 1] -> [L x F]
        W_trf = reshape(W_core, [], F);
    case 'raw'
        % no temporal axis; return empty
        W_trf = [];
end

% --- 6) Package output ---
basis_out = struct();
switch designType
    case 'basis'
        basis_out = struct('B', B, 'lags_ms', lags_ms, ...
                           'centers_ms', centers_ms, 'width_ms', width_ms);
    case 'lagged'
        basis_out = struct('B', [], 'lags_ms', lags_ms, ...
                           'centers_ms', [], 'width_ms', []);
    case 'raw'
        basis_out = struct('B', [], 'lags_ms', [], 'centers_ms', [], 'width_ms', []);
end

out = struct();
out.perf       = rmfield(cvres, 'yhat_oos');
out.cvres      = cvres;
out.yhat_oos   = cvres.yhat_oos;     % aligned to cleaned T only
out.W_trf      = W_trf;              % [] for 'raw'
out.W_basis    = W_basis;            % weights in design space
out.basis      = basis_out;
out.model      = struct('lambda', lam_star, 'mu_X', muX, 'sd_X', sdX, ...
                        'mu_y', muY, 'sd_y', sdY, 'fs', fs, ...
                        'purge_pts', purge_pts, 'inner_metric', opt.InnerMetric, ...
                        'dropped_rows', dropped, 'add_intercept', logical(opt.AddIntercept), ...
                        'DesignType', designType);
out.featNames  = featNames;
out.yhat_final = yhat_final;
end
