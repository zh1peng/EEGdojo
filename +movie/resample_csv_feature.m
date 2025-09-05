function [F_tgt, varNames, t_tgt] = resample_csv_feature(csvPath, varargin)
%resample_csv_feature  Resample a movie-features CSV (with 'seconds') to a target grid.
%
% [F_tgt, varNames, t_tgt] = resample_csv_feature(csvPath, 'Name', Value, ...)
%
% Required:
%   csvPath : path to CSV containing a 'seconds' column + feature columns
%
% Name-Value options:
%   'FsTarget'   : target sampling rate in Hz (required if 'tTarget' not given)
%   'tStart'     : start time (sec). default = min(seconds)
%   'tEnd'       : end time (sec).   default = max(seconds)
%   'tTarget'    : explicit target time vector (sec). Overrides Fs/tStart/tEnd
%   'OffsetSec'  : shift features in time (applied uniformly) (default 0)
%   'HoldVars'   : cellstr of vars to step-hold
%   'InterpVars' : cellstr of vars to continuous-interp
%   'InterpMethod' : 'pchip' (default) or 'linear' for InterpVars
%   'MakePlots'  : true/false sanity plots (default false)
%   'PlotVars'   : subset of vars to plot
%   'Win'        : [t0 t1] seconds to plot (default: first 30s of target)
%   'MaxCols'    : subplot columns for plotting (default 2)
%
% Outputs:
%   F_tgt    : [T_target x F] resampled features
%   varNames : 1xF feature names (order matches columns in F_tgt)
%   t_tgt    : [T_target x 1] target time vector (sec)
%
% % Example (your CSV schema with seconds + emotion counts + detections + AV features)
% csvPath = 'movie_feats.csv';
% % Option A: define just FsTarget; the wrapper builds t_tgt for you
% [F128, names, t128] = resample_csv_feature(csvPath, ...
%     'FsTarget', 128, ...
%     'MakePlots', true, ...
%     'PlotVars', {'spoken_words','faces','brightness','loudness','motion'});

% % Option B: use a custom time grid (e.g., align to EEG epoch times)
% t_eeg = (0:1/200:600)';  % 10 minutes at 200 Hz
% [F200, namesB, t200] = resample_csv_feature(csvPath, ...
%     'tTarget', t_eeg, ...
%     'InterpMethod','pchip', ...
%     'MakePlots', true);

% % Option C: If you want to override how variables are treated:
% holdVars   = {'spoken_words','written_words','letters','has_letters','faces','closeup','body','num_characters', ...
%               'positive','negative','anger','happy','fear','sad','excited'};
% interpVars = {'brightness','saliency','sharpness','vibrance','loudness','tempo','motion'};
% [Fcustom, namesC, tc] = resample_csv_feature(csvPath, ...
%     'FsTarget', 128, ...
%     'HoldVars', holdVars, ...
%     'InterpVars', interpVars, ...
%     'MakePlots', true);


p = inputParser;
p.addParameter('FsTarget', [], @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('tStart', [], @(x)isnumeric(x)&&isscalar(x));
p.addParameter('tEnd',   [], @(x)isnumeric(x)&&isscalar(x));
p.addParameter('tTarget', [], @(x)isnumeric(x)&&isvector(x));
p.addParameter('OffsetSec', 0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('HoldVars', {}, @(x)iscellstr(x)||isempty(x));
p.addParameter('InterpVars', {}, @(x)iscellstr(x)||isempty(x));
p.addParameter('InterpMethod','pchip', @(s)ischar(s)&&ismember(lower(s),{'pchip','linear'}));
p.addParameter('MakePlots', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('PlotVars', {}, @(x)iscellstr(x)||isempty(x));
p.addParameter('Win', [], @(x)isnumeric(x)&&numel(x)==2 || isempty(x));
p.addParameter('MaxCols', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.parse(varargin{:});
opt = p.Results;

% ---------- read CSV ----------
T = readtable(csvPath);
assert(any(strcmpi(T.Properties.VariableNames,'seconds')), 'CSV must have a "seconds" column');
t_src = T.seconds(:);

% Candidate feature names (exclude 'seconds')
varNames = T.Properties.VariableNames;
varNames(strcmpi(varNames,'seconds')) = [];

% ---------- target grid ----------
if ~isempty(opt.tTarget)
    t_tgt = opt.tTarget(:);
else
    assert(~isempty(opt.FsTarget), 'Provide FsTarget or an explicit tTarget.');
    t0 = ~isempty(opt.tStart) * opt.tStart + isempty(opt.tStart) * min(t_src);
    t1 = ~isempty(opt.tEnd)   * opt.tEnd   + isempty(opt.tEnd)   * max(t_src);
    assert(t1 > t0, 'tEnd must be > tStart');
    N  = floor((t1 - t0)*opt.FsTarget) + 1;
    t_tgt = (0:N-1)'/opt.FsTarget + t0;
end
K = numel(t_tgt);

% ---------- decide hold vs interp sets ----------
defaultHold = { ...
    'positive','negative','anger','happy','fear','sad','excited', ... % counts
    'closeup','body','faces','num_characters', ...                    % binary/count
    'spoken_words','written_words','letters','has_letters'            % text counts/flags
    };

% Start from defaults but respect user overrides
if isempty(opt.HoldVars)
    holdVars = intersect(defaultHold, varNames, 'stable');
else
    holdVars = intersect(opt.HoldVars, varNames, 'stable');
end

if isempty(opt.InterpVars)
    interpVars = setdiff(varNames, holdVars, 'stable'); % everything else continuous
else
    interpVars = intersect(opt.InterpVars, varNames, 'stable');
    % Ensure no conflicts
    holdVars   = setdiff(holdVars, interpVars, 'stable');
end

% ---------- resample per variable (via the unit function) ----------
F_tgt = zeros(K, numel(varNames));
% Hold/step variables
for i = 1:numel(holdVars)
    v = holdVars{i};
    F_tgt(:, strcmp(varNames, v)) = resample_one_feature( ...
        t_src, T.(v), t_tgt, 'hold', ...
        'OffsetSec', opt.OffsetSec, ...
        'Fill','ffill');
end
% Continuous variables
for i = 1:numel(interpVars)
    v = interpVars{i};
    F_tgt(:, strcmp(varNames, v)) = resample_one_feature( ...
        t_src, T.(v), t_tgt, 'interp', ...
        'InterpMethod', opt.InterpMethod, ...
        'OffsetSec', opt.OffsetSec, ...
        'Fill','ffill');
end

% ---------- optional sanity plots ----------
if opt.MakePlots
    local_plot(T, t_src + opt.OffsetSec, F_tgt, varNames, t_tgt, opt);
end
end

% ===== local plotting helper =====
function local_plot(Tsrc, t_feat, F, varNames, t_tgt, opt)
    defaultVars = {'spoken_words','faces','written_words','letters', ...
                   'brightness','saliency','loudness','tempo','motion'};
    if isempty(opt.PlotVars)
        vars = intersect(defaultVars, varNames, 'stable');
        if isempty(vars), vars = varNames(1:min(6,numel(varNames))); end
    else
        vars = intersect(opt.PlotVars, varNames, 'stable');
        if isempty(vars), error('Requested PlotVars not found in CSV'); end
    end

    if isempty(opt.Win)
        t0 = t_tgt(1); t1 = min(t_tgt(1)+30, t_tgt(end)); % first 30 s
    else
        t0 = max(opt.Win(1), t_tgt(1));
        t1 = min(opt.Win(2), t_tgt(end));
    end
    idx_tgt = (t_tgt>=t0 & t_tgt<=t1);
    idx_src = (t_feat>=t0 & t_feat<=t1);

    n = numel(vars); ncols = min(opt.MaxCols, n); nrows = ceil(n/ncols);
    figure('Color','w','Name','CSV→Target resample sanity');
    for i = 1:n
        v = vars{i}; ax=subplot(nrows,ncols,i); hold(ax,'on');
        if any(strcmpi(Tsrc.Properties.VariableNames, v))
            plot(ax, t_feat(idx_src), Tsrc.(v)(idx_src), 'o', 'MarkerSize',4, ...
                 'DisplayName','source');
        end
        ci = find(strcmp(varNames, v), 1);
        if ~isempty(ci)
            plot(ax, t_tgt(idx_tgt), F(idx_tgt, ci), '-', 'LineWidth',1.2, ...
                 'DisplayName','resampled');
        end
        grid(ax,'on'); xlim(ax,[t0 t1]);
        title(ax, v, 'Interpreter','none'); xlabel(ax,'Time (s)'); ylabel(ax,'Value');
        legend(ax,'Location','best'); legend(ax,'boxoff');
    end
    sgtitle(sprintf('Resample sanity: %d vars  |  window %.2f–%.2f s', n, t0, t1));
end
