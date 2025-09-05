function y_tgt = resample_one_feature(t_src, y_src, t_tgt, mode, varargin)
%RESAMPLE_ONE_FEATURE  Resample a single feature onto a target time grid.
%
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, mode, 'Name', Value, ...)
%
% Inputs
%   t_src : [T x 1] source time stamps (seconds, irregular OK)
%   y_src : [T x 1] source values (numeric)
%   t_tgt : [K x 1] target time grid (seconds)
%   mode  : 'hold'   -> step-wise/previous-hold (counts/binary/categorical)
%           'interp' -> continuous interpolation (e.g., brightness, loudness)
%
% Name-Value options
%   'InterpMethod' : interpolation for 'interp' mode: 'pchip' (default) | 'linear'
%   'OffsetSec'    : shift source times by this many seconds (default 0)
%   'Fill'         : fill strategy for NaNs: 'ffill' (default) | 'bfill' | 'zero' | 'none'
%
% Visualization (optional)
%   'MakePlot'     : true/false to show a sanity plot (default: false)
%   'Ax'           : axis handle to draw into (default: [])
%   'PlotWindow'   : [t0 t1] seconds to zoom plot (default: [])
%   'SourceStyle'  : struct with fields for the source scatter (e.g., MarkerSize)
%   'TargetStyle'  : struct with fields for the target line   (e.g., LineWidth)
%   'Title'        : custom title string (default auto)
%
% Output
%   y_tgt : [K x 1] resampled values aligned to t_tgt
%
% Notes
% - Handles irregular/duplicate t_src (keeps first occurrence, stable).
% - Robust to NaNs/Infs; final non-finite samples are set to 0 unless Fill='none'.
% - 'hold' produces a boxcar/step function (last value carried forward).
%
% -------------------------------------------------------------------------
% Examples
%
% % Example 1: Simple hold-resample (boxcar counts)
% t_src = (0:0.5:10)';              % source sampled every 0.5 s
% y_src = randi([0 3], size(t_src));% discrete counts
% t_tgt = (0:0.01:10)';             % target grid at 100 Hz
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, 'hold', ...
%     'MakePlot', true);
%
% % Example 2: Continuous interpolation of a smooth feature
% t_src = (0:0.2:5)';               % source sampled at 5 Hz
% y_src = sin(2*pi*0.5*t_src);      % continuous sine wave
% t_tgt = (0:0.01:5)';              % resample to 100 Hz
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, 'interp', ...
%     'InterpMethod','pchip', ...
%     'MakePlot', true, ...
%     'PlotWindow', [0 2]);
%
% % Example 3: Apply offset and zero-fill missing values
% t_src = [0 1 3 4]';               % irregular timestamps
% y_src = [10 NaN 30 40]';          % contains a NaN
% t_tgt = (0:0.5:5)';               % target grid
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, 'interp', ...
%     'OffsetSec', 0.2, ...
%     'Fill','zero', ...
%     'MakePlot', true, ...
%     'Title','Offset + Zero Fill');
%
% % Example 4: Custom plotting styles and axes
% ax = subplot(2,1,1);
% y_tgt = resample_one_feature(t_src, y_src, t_tgt, 'hold', ...
%     'MakePlot', true, 'Ax', ax, ...
%     'SourceStyle', struct('Marker','s','Color','r'), ...
%     'TargetStyle', struct('Color','g','LineWidth',2));
%
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

% Drop NaN times
ok = isfinite(t_src);
t_src = t_src(ok);
y_src = y_src(ok);

% Handle duplicates (keep first occurrence, stable)
[tu, ia] = unique(t_src, 'stable');
if numel(tu) < numel(t_src)
    % Duplicate timestamps found; keep first occurrences.
end
t_src = tu; y_src = y_src(ia);

% If nothing usable, return zeros
if isempty(t_src) || isempty(y_src) || all(~isfinite(y_src))
    y_tgt = zeros(size(t_tgt));
    if opt.MakePlot, local_plot(t_src, y_src, t_tgt, y_tgt, opt, mode); end
    return;
end

% -------------------- resample --------------------
switch lower(mode)
    case 'hold'
        % Step/hold interpolation: carry last value forward
        y = interp1(t_src, y_src, t_tgt, 'previous', 'extrap');
        % Fill any internal gaps
        y = fillmissing(y,'previous');
        y = fillmissing(y,'next');

    case 'interp'
        % Continuous interpolation
        y = interp1(t_src, y_src, t_tgt, opt.InterpMethod, 'extrap');

    otherwise
        error('Unknown mode "%s". Use "hold" or "interp".', mode);
end

% -------------------- post-fill --------------------
switch lower(opt.Fill)
    case 'ffill'
        y = fillmissing(y,'previous');
        y = fillmissing(y,'next');
    case 'bfill'
        y = fillmissing(y,'next');
        y = fillmissing(y,'previous');
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
