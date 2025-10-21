function y_tgt = resample_one_feature(t_src, y_src, t_tgt, mode, varargin)
%RESAMPLE_ONE_FEATURE  Resample a single feature onto a target time grid.
%
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, mode, 'Name', Value, ...)
%
% Inputs
%   t_src : [T x 1] source timestamps in seconds (irregular OK)
%   y_src : [T x 1] source values (numeric/logical)
%   t_tgt : [K x 1] target time grid in seconds (uniform EEG grid, etc.)
%   mode  : 'hold'   -> step-wise (zero-order hold / previous value)
%           'interp' -> continuous interpolation ('pchip' or 'linear')
%
% Name-Value options
%   'InterpMethod' : interpolation for 'interp' mode:
%                    'pchip' (default) | 'linear'
%   'OffsetSec'    : shift source times by this many seconds (default 0)
%   'Fill'         : handle NaNs after resampling:
%                    'ffill' (default) -> forward-fill internal NaNs
%                    'bfill'           -> back-fill internal NaNs
%                    'zero'            -> set NaNs to 0
%                    'none'            -> leave NaNs
%
% Visualization (optional)
%   'MakePlot'     : true/false to show a sanity plot (default: false)
%   'Ax'           : axis handle to draw into (default: [])
%   'PlotWindow'   : [t0 t1] seconds to zoom plot (default: [])
%   'SourceStyle'  : struct for the source scatter (e.g., MarkerSize)
%   'TargetStyle'  : struct for the target line   (e.g., LineWidth)
%   'Title'        : custom title string (default auto)
%
% Output
%   y_tgt : [K x 1] resampled values aligned to t_tgt
%
% Notes
% - Designed for naturalistic annotations:
%     * Impulse-like/binary/categorical → use mode='hold'
%     * Continuous (brightness/loudness/etc.) → mode='interp'
% - By default, this function does NOT extrapolate outside the source time
%   range. For impulses this avoids “turning on” events before they exist.
% - In 'hold' mode we forward-fill *internal* gaps only (no future leakage).
% - Timestamps are sorted; duplicate t_src are collapsed (keeps first).
% - Final non-finite samples are set to 0 unless Fill='none'.
%
% -------------------------------------------------------------------------
% Examples (paste into Command Window)
%
% % Example 1: Simple hold-resample (binary/impulse labels at 10 Hz → 128 Hz)
% t_src = (0:0.1:5)';                      % 10 Hz
% y_src = double(t_src==1.0 | t_src==3.2); % impulse-like
% t_tgt = (0:1/128:5)';                    % 128 Hz EEG grid
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, 'hold', ...
%     'MakePlot', true, 'Title','Impulse → hold (no extrap)');
%
% % Example 2: Continuous interpolation of a smooth feature
% t_src = (0:0.2:5)';                      % 5 Hz
% y_src = sin(2*pi*0.5*t_src);             % continuous
% t_tgt = (0:0.01:5)';                      % 100 Hz
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, 'interp', ...
%     'InterpMethod','pchip', ...
%     'MakePlot', true, 'PlotWindow', [0 2], ...
%     'Title','Continuous → pchip');
%
% % Example 3: Irregular timestamps, offset and zero-fill
% t_src = [0 1 3 4]';  y_src = [10 NaN 30 40]';
% t_tgt = (0:0.5:5)';   % 2 Hz target
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, 'interp', ...
%     'OffsetSec', 0.2, 'Fill','zero', 'MakePlot', true, ...
%     'Title','Offset + zero-fill');
%
% % Example 4: Custom plotting styles and provided axes
% ax = subplot(2,1,1);
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, 'hold', ...
%     'MakePlot', true, 'Ax', ax, ...
%     'SourceStyle', struct('Marker','s','Color','r'), ...
%     'TargetStyle', struct('Color','g','LineWidth',2));
% -------------------------------------------------------------------------

% -------------------- parse options --------------------
p = inputParser;
p.addParameter('InterpMethod','pchip', @(s)ischar(s)&&ismember(lower(s),{'pchip','linear'}));
p.addParameter('OffsetSec',0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Fill','ffill', @(s)ischar(s)&&ismember(lower(s),{'ffill','bfill','zero','none'}));

% viz
p.addParameter('MakePlot',false, @(x)islogical(x)&&isscalar(x));
p.addParameter('Ax',[], @(h) isempty(h) || isgraphics(h,'axes'));
p.addParameter('PlotWindow',[], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
p.addParameter('SourceStyle',struct(), @isstruct);
p.addParameter('TargetStyle',struct(), @isstruct);
p.addParameter('Title','', @(s)ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

% -------------------- normalize/validate --------------------
t_src = t_src(:) + opt.OffsetSec;
y_src = y_src(:);
t_tgt = t_tgt(:);

% Drop NaN/Inf times
ok = isfinite(t_src);
t_src = t_src(ok);
y_src = y_src(ok);

% If nothing usable, return zeros (and optional plot)
if isempty(t_src) || isempty(y_src) || all(~isfinite(y_src))
    y_tgt = zeros(size(t_tgt));
    if opt.MakePlot, local_plot(t_src, y_src, t_tgt, y_tgt, opt, mode); end
    return;
end

% Sort ascending (interp1 expects monotonic x)
[ts, ord] = sort(t_src, 'ascend');
ys = y_src(ord);

% Handle duplicates on the sorted axis (keep first; change to 'last'/'sum' if needed)
[tu, ia] = unique(ts, 'stable');
t_src = tu;
y_src = ys(ia);

% -------------------- resample (no extrap by default) --------------------
% Edge policy: for stick-like features, set values outside source support to zero
edge_left_zero  = true;
edge_right_zero = true;

switch lower(mode)
    case 'hold'
        % Step/hold (previous) WITHOUT extrapolation
        y = interp1(t_src, y_src, t_tgt, 'previous');
        % Internal gaps only: forward-fill (no look-ahead)
        y = fillmissing(y, 'previous');

    case 'interp'
        % Continuous interpolation WITHOUT extrapolation
        y = interp1(t_src, y_src, t_tgt, opt.InterpMethod);

    otherwise
        error('Unknown mode "%s". Use "hold" or "interp".', mode);
end

% -------------------- edges and post-fill --------------------
% Identify target samples outside the source support
left_edge  = t_tgt < t_src(1);
right_edge = t_tgt > t_src(end);

% Apply edge policy (zeros for impulses; change if sustained labels desired)
if edge_left_zero,  y(left_edge)  = 0; end
if edge_right_zero, y(right_edge) = 0; end

% Now address any remaining internal NaNs based on requested Fill policy
switch lower(opt.Fill)
    case 'ffill'   % forward fill internal NaNs only
        y = fillmissing(y, 'previous');
    case 'bfill'   % back fill internal NaNs only
        y = fillmissing(y, 'next');
    case 'zero'
        y(~isfinite(y)) = 0;
    case 'none'
        % leave NaNs as-is
end

% Final cleanup (avoid Infs)
y(~isfinite(y)) = 0;

y_tgt = y;

% -------------------- optional visualization --------------------
if opt.MakePlot
    local_plot(t_src, y_src, t_tgt, y_tgt, opt, mode);
end
end

% ==================== local helpers ====================
function local_plot(t_src, y_src, t_tgt, y_tgt, opt, mode)
    % Axis
    if isempty(opt.Ax)
        fig = figure('Color','w','Name','Resample One Feature'); %#ok<NASGU>
        ax = axes; hold(ax,'on');
    else
        ax = opt.Ax; hold(ax,'on');
    end

    % Styles (defaults + user overrides)
    srcStyle = struct('Marker','o','MarkerSize',5,'LineStyle','none','Color','k','DisplayName','source');
    tgtStyle = struct('LineWidth',1.4,'Color',[0 0.4470 0.7410],'DisplayName','resampled');
    srcStyle = local_merge_struct(srcStyle, opt.SourceStyle);
    tgtStyle = local_merge_struct(tgtStyle, opt.TargetStyle);

    % Plot
    if ~isempty(t_src) && ~isempty(y_src)
        plot(ax, t_src, y_src, srcStyle);
    end
    plot(ax, t_tgt, y_tgt, tgtStyle);

    % Window
    if ~isempty(opt.PlotWindow)
        xlim(ax, [max(opt.PlotWindow(1), min(t_tgt)) , min(opt.PlotWindow(2), max(t_tgt))]);
    else
        xlim(ax, [min(t_tgt) max(t_tgt)]);
    end

    grid(ax,'on');
    xlabel(ax,'Time (s)'); ylabel(ax,'Value');

    ttl = strtrim(string(opt.Title));
    if ttl == ""
        ttl = sprintf('Resample (%s)', lower(mode));
    end
    title(ax, ttl, 'Interpreter','none');

    legend(ax,'Location','best'); legend(ax,'boxoff');
end

function S = local_merge_struct(S, U)
    if isempty(U), return; end
    f = fieldnames(U);
    for i = 1:numel(f)
        S.(f{i}) = U.(f{i});
    end
end
