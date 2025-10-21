function OUT = mtrf_fit(Y, Xraw, fs, varargin)
% MTRF_FIT  TRF-style encoding wrapper around nested_cv_fit with flexible designs.
%
%   OUT = mtrf_fit(Y, Xraw, fs, 'Name', Value, ...)
%
% Inputs
%   Y    : [T x M] response (neural) time series.
%   Xraw : [T x F] predictor(s) (stimulus features).
%   fs   : sampling rate (Hz).
%
% Key Name-Value options (defaults in [brackets])
%   DesignType     : 'raw' | 'lagged' | 'cosine'                      ['lagged']
%   lags_ms        : vector of lags (ms) for 'lagged'/'cosine'        [ -100:10:500 ]
%   centers_ms     : centers for raised-cosine basis (cosine)         []
%   width_ms       : width for raised-cosine basis (cosine)           []
%   zeropad        : logical, for 'lagged' mode only (mTRF-style)     [true]
%   pad            : 'zero' | 'nan' | 'edge' for zeropad=true         ['zero']
%
%   Mode           : 'assess' | 'deploy' | 'both'                     ['assess']
%   Model          : model string or struct for nested_cv_fit         ['ridge']
%   ParamGrid      : hyperparameter grid for nested_cv_fit            [struct('lambda',2.^(0:2:18))]
%   SelectMetric   : 'corr'|'r2'|'mse'|function_handle                ['r2']
%   KOuter         : outer folds (assessment)                         [5]
%   KInner         : inner folds                                      [4]
%   GuardBandSec   : guard band around outer test (seconds)           [0]
%   UseParallel    : true/false                                       [false]
%   Verbose        : 0/1                                              [1]
%   VariableNames  : cellstr/strings for features (optional)          []
%
% Outputs (struct OUT)
%   .design.type        : 'raw'|'lagged'|'cosine'
%   .design.lags_ms     : lag axis used (if applicable)
%   .design.lag_samp    : integer-sample lags (if applicable)
%   .design.B           : basis (cosine) or [] otherwise
%   .design.kept_idx    : rows retained after padding/trim/NaN drop
%
%   .cv                 : full result from nested_cv_fit (assessment/deploy/both)
%
%   .deploy.W_design_z  : weights (design space) learned in z-space (if deploy/both)
%   .deploy.W_design_eff: weights mapped to ORIGINAL units (see Notes)
%   .deploy.TRF_lag_eff : TRF in original units: [L x F x M] for lagged/cosine (if deploy/both)
%   .deploy.axes        : struct with lag_ms, feature_names (if provided)
%
% Notes
%   • Standardization happens *inside* nested_cv_fit with OUTER-TRAIN stats
%     for assessment and ALL-data stats for deploy. No leakage.
%   • ORIGINAL-units effective weights are computed from z-space model as:
%       Y ≈ ((X - muX)./sdX) * Wz .* sdY + muY  ⇒  Weff = diag(1./sdX) * Wz * diag(sdY)
%   • For 'raw' design there is no lag axis; TRF_lag_eff is [].

% OUT = mtrf_fit(Y, X, fs, ...
%   'DesignType','lagged', ...
%   'lags_ms', -100:10:500, ...
%   'zeropad', true, 'pad','edge', ...   % or 'nan' (then rows with NaN are dropped)
%   'Mode','both', ...
%   'Model','ridge', ...
%   'ParamGrid', struct('lambda',2.^(0:2:18)), ...
%   'SelectMetric','corr', ...
%   'KOuter',5, 'KInner',4, 'GuardBandSec', max(abs([-100 500]))/1000, ...
%   'UseParallel', true);

% Raised-cosine basis:
% B = build_cosine(-100:10:500, centers_ms, width_ms);   % optional preview
% OUT = mtrf_fit(Y, X, fs, ...
%   'DesignType','cosine', ...
%   'lags_ms', -100:10:500, 'centers_ms', centers_ms, 'width_ms', width_ms, ...
%   'Mode','both', 'Model','ridge', 'ParamGrid', struct('lambda',2.^(0:2:18)));

% Example 3
% Raw (no temporal expansion):
% OUT = mtrf_fit(Y, X, fs, 'DesignType','raw', 'Mode','assess');


% ---------- Parse options ----------
ip = inputParser;
ip.addParameter('DesignType','lagged',@(s) any(strcmpi(s,{'raw','lagged','cosine'})));
ip.addParameter('lags_ms', -100:10:500, @isnumeric);
ip.addParameter('centers_ms', [], @isnumeric);
ip.addParameter('width_ms', [], @isnumeric);
ip.addParameter('zeropad', true, @islogical);
ip.addParameter('orth', true, @islogical);
ip.addParameter('pad', 'zero', @(s) any(strcmpi(s,{'zero','nan','edge'})));
ip.addParameter('ShowDesign', true, @islogical);
ip.addParameter('Mode','assess');
ip.addParameter('Model','ridge');
ip.addParameter('ParamGrid', struct('lambda', 2.^(0:2:18)));
ip.addParameter('SelectMetric','r2');
ip.addParameter('KOuter',5,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
ip.addParameter('KInner',4,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
ip.addParameter('GuardBandSec',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('UseParallel',false,@islogical);
ip.addParameter('Verbose',1,@(x)isnumeric(x)&&isscalar(x));
ip.addParameter('VariableNames',[],@(v) isstring(v) || iscellstr(v) || isempty(v));
ip.parse(varargin{:});
opt = ip.Results;

[T1, F] = size(Xraw);
[T2, M] = size(Y);
assert(T1==T2, 'mtrf_fit: Xraw and Y must have the same number of rows.');

% ---------- 1) Build design matrix ----------
designType = lower(opt.DesignType);
B = [];
lag_samp = [];
kept_idx = (1:T1).';

switch designType
    case 'raw'
        Xdesign = Xraw;

    case 'lagged'
        [Xlag, kept_idx, lag_samp] = mtrf.lag_design(Xraw, opt.lags_ms, ...
            'units','ms','fs',fs,'zeropad',opt.zeropad,'pad',opt.pad,'bias',false);
        Xdesign = Xlag;
        % If zeropad=false, kept_idx trims edges; if zeropad=true, keep T rows

    case 'cosine'
        % Build raised-cosine basis & convolve features (your functions)
        B = mtrf.build_cosine(opt.lags_ms, opt.centers_ms, opt.width_ms);  % [L x K]
        
        if opt.orth
            [Q,~] = qr(B,0);      % Q: L×K with Q'Q = I
            Xdesign = mtrf.cosine_design(Xraw, Q); 
        else
            Xdesign = mtrf.cosine_design(Xraw, B);
        end

                                     % [T x (F*K)]
        lag_samp = round(opt.lags_ms .* fs ./ 1000);                  % for metadata
        % cosine_design generally returns T rows; handle NaNs below

    otherwise
        error('mtrf_fit: unknown DesignType "%s".', opt.DesignType);
end

if isfield(opt,'ShowDesign')&&opt.ShowDesign
    mtrf.show_design_matrix(Xraw,fs,'DesignType',...
                    designType,'lags_ms',...
                    opt.lags_ms,'zeropad',...
                    opt.zeropad,'pad',opt.pad,...
                    'centers_ms',opt.centers_ms,...
                    'width_ms',opt.width_ms,...
                    'orth',isfield(opt,'orth')&&opt.orth,'FeatureNames',opt.VariableNames,...
                    'Title',sprintf('Design preview (%s)',designType)); 
end


% ---------- 2) Handle NaNs / alignment ----------
% General rule: drop any rows with NaN in Xdesign or Y (works for all designs)
bad = any(isnan(Xdesign),2) | any(isnan(Y),2);
if any(bad)
    Xdesign = Xdesign(~bad, :);
    Y       = Y(~bad, :);
    kept_idx = kept_idx(~bad);
end

% ---------- 3) Call nested_cv_fit ----------
guard_samples = round(opt.GuardBandSec * fs);
cvOUT = mtrf.nested_cv_fit(Y, Xdesign, ...
    'Mode', opt.Mode, ...
    'Model', opt.Model, ...
    'ParamGrid', opt.ParamGrid, ...
    'SelectMetric', opt.SelectMetric, ...
    'KOuter', opt.KOuter, ...
    'KInner', opt.KInner, ...
    'GuardBand', guard_samples, ...
    'UseParallel', opt.UseParallel, ...
    'Verbose', opt.Verbose);

% ---------- 4) Reconstruct TRF (from deploy model only) ----------
W_design_z  = [];
W_design_eff= [];
TRF_lag_eff = [];

if isfield(cvOUT, 'deploy') && ~isempty(cvOUT.deploy)
    % z-space weights (design space)
    if isfield(cvOUT.deploy.model,'W')
        W_design_z = cvOUT.deploy.model.W;  % [(F*K [+1 if intercept]) x M]
    else
        W_design_z = [];
    end

    % Map to ORIGINAL units: Weff = diag(1./sdX)*Wz*diag(sdY)
    px = size(Xdesign, 2);
    if ~isempty(W_design_z) && isfield(cvOUT.deploy,'preproc')
        muX = cvOUT.deploy.preproc.muX;
        sdX = cvOUT.deploy.preproc.sdX;
        muY = cvOUT.deploy.preproc.muY; %#ok<NASGU> (kept for clarity)
        sdY = cvOUT.deploy.preproc.sdY;

        % If we added an intercept, its scaling is: (1/sdX_bias)*Wz * sdY,
        % but our bias col has sdX≈0 ⇒ we set sdX_bias=1 in nested_cv_fit.
        scaleX = 1 ./ sdX(:);    % [px x 1]
        scaleY = sdY(:).';       % [1 x M]
        W_design_eff = (W_design_z .* scaleX) .* scaleY;  % broadcasting
    end

    % Rebuild lag-domain TRF (original units) if we can
    if ~isempty(W_design_eff)
        % Remove intercept before reshaping, if present
            Wcore = W_design_eff;

        switch designType
            case 'lagged'
                L = numel(opt.lags_ms);
                Fused = size(Wcore,1) / L;
                if mod(size(Wcore,1), L) ~= 0
                    warning('mtrf_fit: shape mismatch for lagged reshape; skipping TRF.');
                else
                    % Columns are lag-major: [lag1 x F, lag2 x F, ...]
                    % Reshape as [F, L, M] then permute -> [L, F, M]
                    WFLM = reshape(Wcore, [Fused, L, M]);   % [F x L x M]
                    TRF_lag_eff = permute(WFLM, [2, 1, 3]); % [L x F x M]
                end

            case 'cosine'
                % Wcore is [(F*K) x M]. Split into K per feature, then multiply by B (LxK)
                L = size(B,1); K = size(B,2);
                Fused = size(Wcore,1) / K;
                if mod(size(Wcore,1), K) ~= 0
                    warning('mtrf_fit: shape mismatch for cosine reshape; skipping TRF.');
                else
                    TRF_lag_eff = zeros(L, Fused, M, 'like', Wcore);
                    for m = 1:M
                        WfK = reshape(Wcore(:,m), K, Fused);   % [K x F]
                        TRF_lag_eff(:,:,m) = B * WfK;          % [L x F]
                    end
                end

            otherwise
                % raw: no lag axis
                TRF_lag_eff = [];
        end
    end
end

% ---------- 5) Package output ----------
OUT = struct();
OUT.design = struct( ...
    'type', designType, ...
    'lags_ms', iff(any(strcmpi(designType,{'lagged','cosine'})), opt.lags_ms, []), ...
    'lag_samp', lag_samp, ...
    'B', iff(strcmpi(designType,'cosine'), B, []), ...
    'kept_idx', kept_idx, ...
    'fs', fs ...
);

if ~isempty(opt.VariableNames)
    OUT.design.feature_names = cellstr(opt.VariableNames);
else
    OUT.design.feature_names = arrayfun(@(i) sprintf('feat%02d',i), 1:F, 'uni',0);
end

OUT.cv = cvOUT;

if ~isempty(W_design_z)
    OUT.deploy = struct();
    OUT.deploy.W_design_z   = W_design_z;
    OUT.deploy.W_design_eff = W_design_eff;
    OUT.deploy.TRF_lag_eff  = TRF_lag_eff;
    OUT.deploy.axes = struct();
    if any(strcmpi(designType,{'lagged','cosine'}))
        OUT.deploy.axes.lag_ms = opt.lags_ms(:);
    else
        OUT.deploy.axes.lag_ms = [];
    end
    OUT.deploy.axes.feature_names = OUT.design.feature_names;
else
    OUT.deploy = struct(); % empty if Mode='assess' only
end

end

% ---------------------------- small utils ----------------------------
function v = iff(cond, a, b)
if cond, v = a; else, v = b; end
end

