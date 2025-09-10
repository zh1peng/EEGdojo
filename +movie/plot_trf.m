function H = plot_trf(TRFresList, varargin)
% PLOT_TRF  Plot TRFs (mean±SE or individual) from 1+ multivariate group_mtrf_fit results.
%
%   H = plot_trf(TRFresList, 'predictors', {...}, 'from','x'|'cov'|'all', ...)
%
% INPUT
%   TRFresList : One struct OR a cell array / vector of structs returned by group_mtrf_fit.
%                Each struct must contain:
%                  .t (1 x nLags)
%                  .predNames (1 x nPred cell)
%                  .xIdx, .covIdx
%                  .W  (S x nPred x nLags)
%                  .mu (nLags x nPred)
%                  .se (nLags x nPred)   (optional)
%
% KEY OPTIONS
%   'predictors'     : predictor selection (names or indices). If empty, uses 'from'.
%   'from'           : default selection when 'predictors' is empty: 'x'|'cov'|'all'  ['x']
%   'modelLabels'    : cellstr labels for models in legend (default {'M1','M2',...})
%   'plotMode'       : 'mean' (default) or 'individual'
%   'showSE'         : true/false (default true; mean mode only)
%   'timeWindow'     : [tmin tmax] ms (default: overlap across models)
%   'axes'           : target axes (default: new figure)
%   'alpha'          : 0..1 face alpha for SE band (mean) or line lightening (individual)
%   'lineWidth'      : mean line width (default 2)
%   'indivLineWidth' : individual line width (default 0.75)
%   'lineStyles'     : cycles across series (default {'-'})
%   'colors'         : Kx3 colormap to cycle per series (default: axes ColorOrder)
%   'seriesProps'    : per-series graphics overrides (cell of NV cells, or struct array)
%   'legend'         : 'on'|'off' (default 'on')
%   'legendLocation' : legend location (default 'best')
%   'xlabel','ylabel','title' : labels (defaults provided)
%
% OUTPUT
%   H.ax, H.meanLines, H.bands, H.indiv, H.legend
%
% NOTE
%   A "series" = one (model, predictor) pair. This function expands
%   your inputs into a list of series and plots them together.

% EXAMPLE
% TRFres=movie.group_mtrf_fit(Ys, X, cov, xNames, covNames,250);
% H = movie.plot_trf(TRFres, 'from','x', 'plotMode','mean', 'timeWindow',[0 500]);
% H = movie.plot_trf(TRFres, 'predictors', {'brightness','motion'});

% TRFres1=movie.group_mtrf_fit(Ys, F_tgt(1:end-1, 1), cov, varNames(1), covNames,250);
% TRFres2=movie.group_mtrf_fit(Ys, F_tgt(1:end-1, 2), cov, varNames(2), covNames,250);
% H = movie.plot_trf({TRFres1, TRFres2}, ...
%              'predictors', {'positive','negative'}, ...
%              'modelLabels', {'x1+cov','x2+cov'}, ...
%              'plotMode','mean', 'showSE', true, 'timeWindow',[0 500]);

% ---------- Parse options ----------
ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('predictors', [], @(x)isnumeric(x) || isstring(x) || iscellstr(x));
ip.addParameter('from', 'x', @(s)any(strcmpi(s,{'x','cov','all'})));
ip.addParameter('modelLabels', [], @(c)isempty(c) || iscellstr(c) || isstring(c));

ip.addParameter('plotMode', 'mean', @(s)any(strcmpi(s,{'mean','individual'})));
ip.addParameter('showSE', true, @islogical);
ip.addParameter('timeWindow', [], @(v)isnumeric(v) && (isempty(v) || numel(v)==2));
ip.addParameter('axes', [], @(h)isempty(h) || isgraphics(h,'axes'));
ip.addParameter('alpha', [], @(x)isempty(x) || (isscalar(x) && x>=0 && x<=1));
ip.addParameter('lineWidth', 2, @(x)isscalar(x) && x>0);
ip.addParameter('indivLineWidth', 0.75, @(x)isscalar(x) && x>0);
ip.addParameter('lineStyles', {'-'}, @(c)iscellstr(c) && ~isempty(c));
ip.addParameter('colors', [], @(x)isempty(x) || (isnumeric(x) && size(x,2)==3));
ip.addParameter('seriesProps', [], @(x)isempty(x) || iscell(x) || isstruct(x));
ip.addParameter('legend', 'on', @(s)any(strcmpi(s,{'on','off'})));
ip.addParameter('legendLocation', 'best', @(s)ischar(s) || isstring(s));
ip.addParameter('xlabel', 'Time lag (ms)', @(s)ischar(s) || isstring(s));
ip.addParameter('ylabel', 'TRF weight', @(s)ischar(s) || isstring(s));
ip.addParameter('title', '', @(s)ischar(s) || isstring(s));
ip.parse(varargin{:});
opt = ip.Results;

% ---------- Normalize list of models ----------
if isstruct(TRFresList) && ~isscalar(TRFresList)
    models = num2cell(TRFresList);
elseif isstruct(TRFresList)
    models = {TRFresList};
elseif iscell(TRFresList)
    models = TRFresList;
else
    error('TRFresList must be a struct or a cell array of structs.');
end
nM = numel(models);
if isempty(opt.modelLabels)
    modelLabels = arrayfun(@(i)sprintf('M%d',i), 1:nM, 'uni',0);
else
    modelLabels = cellstr(opt.modelLabels);
    assert(numel(modelLabels)==nM, 'modelLabels length must match #models.');
end

% ---------- Validate required fields & gather names ----------
nameLists = cell(1,nM);
for m = 1:nM
    mustHave = {'t','predNames','xIdx','covIdx','W','mu'};
    for k = 1:numel(mustHave)
        assert(isfield(models{m}, mustHave{k}), 'Missing field "%s" in model %d.', mustHave{k}, m);
    end
    nameLists{m} = string(models{m}.predNames(:)).';  % store as string row
end

% ---------- Resolve predictor names to plot ----------
if isempty(opt.predictors)
    wantNames = strings(0);
    for m = 1:nM
        switch lower(opt.from)
            case 'x',   idx = models{m}.xIdx(:).';
            case 'cov', idx = models{m}.covIdx(:).';
            case 'all', idx = 1:numel(nameLists{m});
        end
        wantNames = union(wantNames, nameLists{m}(idx), 'stable');
    end
else
    if isnumeric(opt.predictors)
        idx = opt.predictors(:).';
        nPred1 = numel(nameLists{1});
        if any(idx<1 | idx>nPred1)
            error('predictors indices out of range 1..%d (refers to 1st model).', nPred1);
        end
        wantNames = nameLists{1}(idx);
    else
        wantNames = string(cellstr(opt.predictors));
    end
end
if isempty(wantNames)
    error('No predictors selected to plot.');
end

% ---------- Expand into plot series (model × predictor) ----------
series = struct('name',[],'label',[],'t',[],'mu',[],'se',[],'W',[]);
series(1) = []; % start empty
for m = 1:nM
    t  = models{m}.t(:).';
    nm_strings = nameLists{m};       % string row
    nm_cell    = cellstr(nm_strings);% cellstr row, robust for strcmpi
    for k = 1:numel(wantNames)
        wn =  wantNames(k);
        wn_char = char(wn); % robust target for strcmpi
        j = find(strcmpi(nm_strings, wn_char), 1, 'first');  % <-- FIXED
        if isempty(j), continue; end  % predictor not in this model
        % dims in your TRFres: mu/se = [nLags x nPred], W = [S x nPred x nLags]
        mu_m = models{m}.mu(:, j);   % [nLags x 1]
        se_m = [];
        if isfield(models{m},'se') && ~isempty(models{m}.se)
            se_m = models{m}.se(:, j); % [nLags x 1]
        end
        W_m  = squeeze(models{m}.W(:, j, :));  % [S x nLags]
        series(end+1) = struct( ... %#ok<AGROW>
            'name',  string(wn), ...
            'label', sprintf('%s: %s', modelLabels{m}, wn_char), ...
            't',     t, ...
            'mu',    mu_m(:).', ...
            'se',    se_m(:).', ...
            'W',     W_m );
    end
end
K = numel(series);
if K==0
    error('None of the requested predictors were found in the provided models.');
end

% ---------- Axes, colors, styles ----------
if isempty(opt.axes) || ~isvalid(opt.axes)
    figure('Color','w'); ax = axes; hold(ax,'on'); box(ax,'on');
else
    ax = opt.axes; hold(ax,'on'); box(ax,'on');
end

% Choose colormap properly
if isempty(opt.colors)
    C = get(ax,'ColorOrder');   % N×3
else
    C = opt.colors;             % K×3
end
LS = opt.lineStyles;

% Alpha defaults
if strcmpi(opt.plotMode,'mean')
    if isempty(opt.alpha), alphaBand = 0.15; else, alphaBand = opt.alpha; end
else
    if isempty(opt.alpha), alphaInd = 0.35; else, alphaInd = opt.alpha; end %#ok<NASGU>
end

% ---------- Time window ----------
if isempty(opt.timeWindow)
    tmins = arrayfun(@(s)min(s.t), series);
    tmaxs = arrayfun(@(s)max(s.t), series);
    tWin = [max(tmins) min(tmaxs)];
else
    tWin = opt.timeWindow(:).';
end
if tWin(1) >= tWin(2), tWin = []; end

% ---------- Prepare outputs ----------
H.ax        = ax;
H.meanLines = gobjects(1,K);
H.bands     = gobjects(1,K);
H.indiv     = cell(1,K);
H.legend    = gobjects(1,1);

% ---------- Plot ----------
warned = false;
for i = 1:K
    s = series(i);
    t = s.t(:).';
    mask = true(size(t));
    if ~isempty(tWin), mask = (t>=tWin(1)) & (t<=tWin(2)); end
    if ~any(mask)
        warning('No samples in window for series %d (%s).', i, s.label);
        continue;
    end

    col = C(mod(i-1, size(C,1))+1, :);
    ls  = LS{mod(i-1, numel(LS))+1};
    extraNV = extract_series_props(opt.seriesProps, i);

    switch lower(opt.plotMode)
        case 'mean'
            mu = s.mu;
            if numel(mu) ~= numel(t) && ~warned
                warning('Length mismatch (t vs mu); truncating to window.'); warned = true;
            end
            mu = mu(mask); tx = t(mask);

            % SE shading (if available)
            if opt.showSE && ~isempty(s.se)
                se = s.se(mask);
                H.bands(i) = fill(ax, [tx, fliplr(tx)], [mu+se, fliplr(mu-se)], col, ...
                    'FaceAlpha', alphaBand, 'EdgeColor','none', 'HandleVisibility','off');
            else
                H.bands(i) = gobjects(1);
            end

            H.meanLines(i) = plot(ax, tx, mu, 'Color', col, 'LineStyle', ls, ...
                'LineWidth', opt.lineWidth, 'DisplayName', s.label);

            if ~isempty(extraNV)
                set(H.meanLines(i), extraNV{:});
                if ishghandle(H.bands(i))
                    k = find(strcmpi(extraNV(1:2:end),'Color'),1);
                    if ~isempty(k), set(H.bands(i),'FaceColor', extraNV{2*k}); end
                end
            end

        case 'individual'
            W = s.W; % [S x nLags]
            if size(W,2) ~= numel(t) && ~warned
                warning('Length mismatch (t vs W); truncating to window.'); warned = true;
            end
            Wm = W(:, mask); tx = t(mask);
            nSubs = size(Wm,1);
            H.indiv{i} = gobjects(1,nSubs);
            % emulate "alpha" by lightening color (lines do not support RGBA)
            liteCol = lighten_color(col, iff(isempty(opt.alpha), 0.6, opt.alpha));
            for ss = 1:nSubs
                H.indiv{i}(ss) = plot(ax, tx, Wm(ss,:), 'Color', liteCol, ...
                    'LineStyle', ls, 'LineWidth', opt.indivLineWidth, ...
                    'HandleVisibility','off');
                if ~isempty(extraNV), set(H.indiv{i}(ss), extraNV{:}); end
            end
            % invisible mean for legend
            H.meanLines(i) = plot(ax, tx, mean(Wm,1,'omitnan'), 'Color', col, ...
                'LineStyle', ls, 'LineWidth', max(1,opt.lineWidth-0.5), ...
                'DisplayName', s.label, 'Visible','off');
            H.bands(i) = gobjects(1);
    end
end

% ---------- Labels, title, legend ----------
xlabel(ax, opt.xlabel);
ylabel(ax, opt.ylabel);
if strlength(string(opt.title))>0
    ttl = string(opt.title);
else
    ttl = 'TRF (mean \pm SE) across subjects';
end
title(ax, ttl);

if strcmpi(opt.legend,'on')
    hasName = arrayfun(@(h) ishghandle(h) && isprop(h,'DisplayName') && ~strcmp(get(h,'Visible'),'off'), H.meanLines);
    L = H.meanLines(hasName);
    if ~isempty(L)
        H.legend = legend(ax, L, 'Location', opt.legendLocation, 'Box','off', 'Interpreter','none');
    end
end
end

% ================= helpers =================
function extraNV = extract_series_props(seriesProps, i)
extraNV = {};
if isempty(seriesProps), return; end
if iscell(seriesProps)
    if numel(seriesProps) >= i && ~isempty(seriesProps{i}), extraNV = seriesProps{i}; end
elseif isstruct(seriesProps)
    if numel(seriesProps) >= i
        f = fieldnames(seriesProps(i));
        extraNV = cell(1, 2*numel(f));
        for k=1:numel(f)
            extraNV{2*k-1} = f{k}; extraNV{2*k} = seriesProps(i).(f{k});
        end
    end
end
end

function out = iff(c,a,b)
if c, out = a; else, out = b; end
end

function c2 = lighten_color(c, amount)
% LIGHTEN_COLOR  Return a lighter color by mixing with white.
amount = max(0,min(1,amount));
c2 = 1 - (1 - c) * (1 - amount);
end
