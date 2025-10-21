function [Xdesign, meta] = show_design_matrix(Xraw, fs, varargin)
% SHOW_DESIGN_MATRIX  Build and visualize TRF design matrices.
%
%   [Xdesign, meta] = show_design_matrix(Xraw, fs, 'Name', Value, ...)
%
% Inputs
%   Xraw : [T x F] stimulus/features (rows=time, cols=features)
%   fs   : sampling rate (Hz)
%
% Name-Value options (defaults in brackets)
%   'DesignType'   : 'raw' | 'lagged' | 'cosine'                 ['lagged']
%   'lags_ms'      : vector of lags (ms)                          [-100:10:500]
%   'zeropad'      : logical (lagged)                             [true]
%   'pad'          : 'zero'|'nan'|'edge' (lagged, when zeropad)   ['zero']
%   'centers_ms'   : cosine centers (ms)                          []
%   'width_ms'     : cosine width (ms)                            []
%   'orth'         : logical, orthonormalize cosine basis via QR  [false]
%   'bias'         : prepend bias col (if your builders support)  [false]
%   'FeatureNames' : cellstr/strings for labels                   []
%   'MaxRowsPreview': #rows to preview in images                  [1500]
%   'MaxColsCorr'  : #columns for corr heatmap                    [60]
%   'ClipPrct'     : percentile for symmetric color clipping      [98]
%   'Title'        : figure title                                 ['Design preview']
%
% Outputs
%   Xdesign : [T' x D] design matrix (after padding/trim; may equal T)
%   meta    : struct with details (see fields below)

% -------------------- parse --------------------
ip = inputParser;
ip.addParameter('DesignType','lagged',@(s) any(strcmpi(s,{'raw','lagged','cosine'})));
ip.addParameter('lags_ms', -100:10:500, @isnumeric);
ip.addParameter('zeropad', true, @islogical);
ip.addParameter('pad', 'zero', @(s) any(strcmpi(s,{'zero','nan','edge'})));
ip.addParameter('centers_ms', [], @isnumeric);
ip.addParameter('width_ms', [], @isnumeric);
ip.addParameter('orth', false, @islogical);
ip.addParameter('bias', false, @islogical);
ip.addParameter('FeatureNames', [], @(v) iscellstr(v) || isstring(v) || isempty(v));
ip.addParameter('MaxRowsPreview', 1500, @(x)isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('MaxColsCorr', 60, @(x)isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('ClipPrct', 98, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=100);
ip.addParameter('Title','Design preview',@(s)ischar(s)||isstring(s));
ip.parse(varargin{:});
opt = ip.Results;

[T, F] = size(Xraw);
designType = lower(opt.DesignType);
B = []; kept_idx = (1:T).'; lag_samp = [];

% -------------------- build design --------------------
switch designType
    case 'raw'
        Xdesign = Xraw;

    case 'lagged'
        [Xlag, kept_idx, lag_samp] = mtrf.lag_design(Xraw, opt.lags_ms, ...
            'units','ms','fs',fs,'zeropad',opt.zeropad,'pad',opt.pad,'bias',false);
        Xdesign = Xlag;

    case 'cosine'
        B = mtrf.build_cosine(opt.lags_ms, opt.centers_ms, opt.width_ms);  % [L x K]
        if opt.orth
            [Q,~] = qr(B,0);  % orthonormal columns
            B = Q;
        end
        Xdesign = mtrf.cosine_design(Xraw, B);                             % [T x (F*K)]
        lag_samp = round(opt.lags_ms .* fs ./ 1000);

    otherwise
        error('show_design_matrix: unknown DesignType "%s".', opt.DesignType);
end

[Td, D] = size(Xdesign);

% -------------------- metadata & column map --------------------
meta = struct();
meta.type        = designType;
meta.lags_ms     = iff(any(strcmpi(designType,{'lagged','cosine'})), opt.lags_ms(:), []);
meta.lag_samp    = iff(any(strcmpi(designType,{'lagged','cosine'})), lag_samp, []);
meta.kept_idx    = kept_idx;
meta.B           = iff(strcmpi(designType,'cosine'), B, []);
meta.orthonormal = iff(strcmpi(designType,'cosine'), is_orthonormal(B), false);
meta.F           = F; meta.T = Td; meta.D = D;

if strcmpi(designType,'lagged')
    L = numel(opt.lags_ms);
    colmap = zeros(L, F);
    for li=1:L, for f=1:F, colmap(li,f) = (li-1)*F + f; end, end
    meta.colmap.lag_f = colmap;  % lag-major
    meta.order_note = 'lag-major: [lag1 × F, lag2 × F, ...]';
elseif strcmpi(designType,'cosine')
    K = size(B,2);
    colmap = zeros(K, F);
    for k=1:K, for f=1:F, colmap(k,f) = (k-1)*F + f; end, end
    meta.colmap.basis_f = colmap; % basis-major
    meta.order_note = 'basis-major: [b1 × F, b2 × F, ...]';
else
    meta.colmap = [];
    meta.order_note = 'raw: columns = features';
end

if isempty(opt.FeatureNames)
    meta.feature_names = arrayfun(@(i)sprintf('feat%02d',i), 1:F, 'uni',0);
else
    meta.feature_names = cellstr(opt.FeatureNames);
end

% -------------------- visualization --------------------
fig = figure('Name', char(opt.Title), 'Color','w');
tiledlayout(fig, 2, 2, 'Padding','compact','TileSpacing','compact');

% panel 1: preview image of design (uniform sampling over time)
nexttile(1);
rows = min(opt.MaxRowsPreview, Td);
rows_idx = unique(round(linspace(1, Td, rows)));
Xprev = Xdesign(rows_idx, :);

% robust symmetric color limits around 0 (no fixed [-1,1])
c = prctile(abs(Xprev(:)), opt.ClipPrct);
if c <= 0 || ~isfinite(c)
    imagesc(Xprev);
else
    imagesc(Xprev, [-c c]);
end
axis tight; colorbar; colormap(parula);
xlabel('Design columns'); ylabel('Time (sampled rows)');
title(sprintf('%s design (uniform preview %d/%d rows)', upperFirst(designType), numel(rows_idx), Td));

% draw separators for groups (lags or bases)
hold on;
if strcmpi(designType,'lagged')
    L = numel(opt.lags_ms);
    for li = 1:L-1
        xline(li*F+0.5, ':', 'Color',[0 0 0 0.25]);
    end
    xticks( (0.5 + F*(0:L-1)) + F/2 );
    xticklabels(string(opt.lags_ms));
    xlabel('Lag (ms) groups');
elseif strcmpi(designType,'cosine')
    K = size(B,2);
    for k = 1:K-1
        xline(k*F+0.5, ':', 'Color',[0 0 0 0.25]);
    end
    xticks( (0.5 + F*(0:K-1)) + F/2 );
    xticklabels("b"+string(1:K));
    xlabel('Basis groups');
end
hold off;

% panel 2: column correlation (subset) — no forced [-1,1] limits
nexttile(2);
cols_corr = min(opt.MaxColsCorr, D);
C = corr(Xdesign(:,1:cols_corr), Xdesign(:,1:cols_corr), 'Rows','pairwise');
imagesc(C); axis image; colormap(redblue()); colorbar;
title(sprintf('Column correlation (first %d/%d cols)', cols_corr, D));
xlabel('Columns'); ylabel('Columns');

% panel 3: Gram of basis (cosine only) to see orthonormality
nexttile(3);
if strcmpi(designType,'cosine')
    G = B.'*B;
    imagesc(G); axis image; colorbar; colormap(parula);
    title(sprintf('B^T B  (orthonormal=%d)', meta.orthonormal));
    xlabel('Basis'); ylabel('Basis');
else
    axis off; text(0.1,0.5,'(No basis for RAW/LAGGED)','FontAngle','italic');
end

% panel 4: simple spy/sparsity
nexttile(4);
spy(abs(Xdesign(rows_idx,:))>0); axis tight;
title('Sparsity (preview)'); xlabel('Columns'); ylabel('Rows');

sgtitle(sprintf('%s — %s | %d rows × %d cols | %s', ...
    char(opt.Title), upperFirst(designType), Td, D, meta.order_note));

end

% ===================== helpers =====================
function tf = is_orthonormal(B)
tf = false;
if isempty(B), return; end
E = B.'*B - eye(size(B,2));
tf = norm(E,'fro') < 1e-8;
end

function out = iff(c,a,b)
if c, out=a; else, out=b; end
end

function s = upperFirst(s)
s = char(s); if isempty(s), return; end
s(1) = upper(s(1));
end

function cmap = redblue()
% simple diverging colormap (blue-white-red)
n = 256; r = (0:n-1)'/max(n-1,1);
cmap = [r, zeros(n,1), flipud(r)];
end
