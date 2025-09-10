function [H, Stats] = plot_trf_diff(TRFin, varargin)
% PLOT_TRF_CONTRAST2  Contrast two TRFs (within a model or across models) with cluster/TFCE stats.
%
%   [H, Stats] = plot_trf_contrast2(TRFin, Name, Value, ...)
%
% Two modes:
%   1) 'within' (default): TRFin is ONE TRFres struct. You must specify:
%         'pos' : predictor selector (name or index) for the POSITIVE condition
%         'neg' : predictor selector (name or index) for the NEGATIVE condition
%      We will extract those two predictors from the SAME TRFres.
%
%   2) 'across': TRFin is a 1x2 cell {TRFres_pos, TRFres_neg} (or 1x2 struct array).
%         'posSel' (optional): predictor selector (name or index) from TRFres_pos (default: 1)
%         'negSel' (optional): predictor selector (name or index) from TRFres_neg (default: 1)
%      We will extract one predictor from EACH TRFres.
%
% INPUT TRFres fields (as in your plot_trf contract):
%   .t (1 x nLags)
%   .predNames (1 x nPred cell)
%   .W  (S x nPred x nLags)    -- subject TRFs per predictor
%   .mu (nLags x nPred)        -- optional (display only)
%   .se (nLags x nPred)        -- optional (display only)
%
% KEY OPTIONS
%   'mode'           : 'within' (default) | 'across'
%   'pos','neg'      : selector for within-mode (name or index)
%   'posSel','negSel': selector for across-mode (defaults: 1, 1)
%   'timeWindow'     : [tmin tmax] ms (default: full overlap)
%   'alpha'          : cluster-forming alpha (default 0.05)
%   'nPerm'          : permutations (default 1000)
%   'tail'           : 0 two-sided | +1 pos | -1 neg (default 0)
%   'minClusterPts'  : min contiguous samples per cluster (default 1)
%   'method'         : 'cluster' (default) | 'tfce'
%   'tfceParams'     : struct with fields H[2], E[0.5], dh[0.1], minRun[1]
%   'showOriginal'   : true/false (default true) show mean±SE of each
%   'diffColor'      : [1x3] RGB (default [0.30 0.30 0.90])
%   'posColor'       : [1x3] RGB (default [0.20 0.70 0.20])
%   'negColor'       : [1x3] RGB (default [0.80 0.20 0.20])
%   'highlightColor' : [1x3] RGB (default auto-lighten from diffColor)
%   'bandAlpha'      : SE fill alpha (default 0.15)
%   'lineWidth'      : line width (default 2)
%   'axes'           : target axes (default new figure)
%   'legend'         : 'on'|'off' (default 'on')
%   'xlabel'         : default 'Time lag (ms)'
%   'ylabel'         : default '\Delta TRF (pos - neg)'
%   'title'          : default 'Difference (mean \pm SE)'
%   'drawZeroLine'   : true/false (default true)
%
% OUTPUTS
%   H     : handles (.ax, .origBands, .origLines, .diffBand, .diffLine, .zero, .sig)
%   Stats : stats struct (cluster or TFCE path; includes .p_corr and masks)
%
% EXAMPLES
%   % A) Same-model contrast
%    Compare 'positive' vs 'negative' inside one TRFres:
% [H,S] = plot_trf_diff(TRFres, ...
%     'mode','within', ...
%     'pos','positive', 'neg','negative', ...
%     'timeWindow',[0 500], 'method','cluster', 'nPerm', 2000);

% B) Across-model contrast
% %    Compare predictor #1 from TRFres1 vs predictor 'negative' from TRFres2:
% [H,S] = plot_trf_diff({TRFres1, TRFres2}, ...
%     'mode','across', ...
%     'posSel',1, 'negSel','negative', ...
%     'timeWindow',[0 500], 'method','tfce', 'nPerm', 3000);

% ---------------- Parse options ----------------
P = inputParser;
P.addParameter('mode','within',@(s) any(strcmpi(s,{'within','across'})));

% selectors
P.addParameter('pos',[],@(x) isnumeric(x) || ischar(x) || isstring(x));
P.addParameter('neg',[],@(x) isnumeric(x) || ischar(x) || isstring(x));
P.addParameter('posSel',1,@(x) isnumeric(x) || ischar(x) || isstring(x));
P.addParameter('negSel',1,@(x) isnumeric(x) || ischar(x) || isstring(x));

% stats/plot
P.addParameter('timeWindow',[],@(v) isnumeric(v) && (isempty(v) || numel(v)==2));
P.addParameter('alpha',0.05,@(x) isscalar(x) && x>0 && x<1);
P.addParameter('nPerm',1000,@(x) isscalar(x) && x>=10);
P.addParameter('tail',0,@(x) isscalar(x) && any(x==[0,1,-1]));
P.addParameter('minClusterPts',1,@(x) isscalar(x) && x>=1);
P.addParameter('method','cluster',@(s) ischar(s) || isstring(s));
P.addParameter('tfceParams',struct(),@isstruct);
P.addParameter('showOriginal',true,@islogical);

% colors/styling
P.addParameter('diffColor',[0.30 0.30 0.90],@(x) isnumeric(x) && isequal(size(x),[1,3]));
P.addParameter('posColor',[0.20 0.70 0.20],@(x) isnumeric(x) && isequal(size(x),[1,3]));
P.addParameter('negColor',[0.80 0.20 0.20],@(x) isnumeric(x) && isequal(size(x),[1,3]));
P.addParameter('highlightColor',[],@(x) isempty(x) || (isnumeric(x) && isequal(size(x),[1,3])));
P.addParameter('bandAlpha',0.15,@(x) isscalar(x) && x>=0 && x<=1);
P.addParameter('lineWidth',2,@(x) isscalar(x) && x>0);

% axes/labels
P.addParameter('axes',[],@(x) isempty(x) || isgraphics(x,'axes'));
P.addParameter('legend','on',@(s) any(strcmpi(s,{'on','off'})));
P.addParameter('xlabel','Time lag (ms)',@(s) ischar(s) || isstring(s));
P.addParameter('ylabel','\Delta TRF (pos - neg)',@(s) ischar(s) || isstring(s));
P.addParameter('title','Difference (mean \pm SE)',@(s) ischar(s) || isstring(s));
P.addParameter('drawZeroLine',true,@islogical);

P.parse(varargin{:});
opt = P.Results;
opt.method = lower(string(opt.method));

if isempty(opt.highlightColor)
    opt.highlightColor = min(opt.diffColor*1.5, 1); % auto-lighten
end

% ---------------- Extract two TRFs (pos,neg) into [S x T] matrices ----------------
switch lower(opt.mode)
    case 'within'
        assert(isstruct(TRFin) && isscalar(TRFin), ...
            'Within-mode requires a SINGLE TRFres struct as first arg.');
        assert(~isempty(opt.pos) && ~isempty(opt.neg), ...
            'Within-mode: you must provide ''pos'' and ''neg'' selectors.');
        TRFpos = select_trf(TRFin, opt.pos);
        TRFneg = select_trf(TRFin, opt.neg);

    case 'across'
        assert(iscell(TRFin) || (isstruct(TRFin) && numel(TRFin)==2), ...
            'Across-mode requires {TRFres_pos, TRFres_neg} or 1x2 struct array.');
        if iscell(TRFin), A = TRFin{1}; B = TRFin{2}; else, A = TRFin(1); B = TRFin(2); end
        TRFpos = select_trf(A, opt.posSel);
        TRFneg = select_trf(B, opt.negSel);

    otherwise
        error('Unknown mode: %s', opt.mode);
end

% ---------------- Basic checks & window ----------------
t1 = TRFpos.t(:).'; t2 = TRFneg.t(:).';
assert(numel(t1)==numel(t2) && all(t1==t2), 't-axes differ between pos/neg.');
t = t1;

assert(size(TRFpos.W,2)==numel(t) && size(TRFneg.W,2)==numel(t), 'W/t size mismatch.');
assert(size(TRFpos.W,1)==size(TRFneg.W,1), 'Different #subjects in pos vs. neg.');

if isempty(opt.timeWindow)
    maskWin = true(size(t));
else
    maskWin = (t>=opt.timeWindow(1)) & (t<=opt.timeWindow(2));
    assert(any(maskWin), 'Requested timeWindow has no samples.');
end
tx = t(maskWin);

% ---------------- Build difference & pointwise stats ----------------
N   = size(TRFpos.W,1);
D   = TRFpos.W(:,maskWin) - TRFneg.W(:,maskWin); % N x T
muD = mean(D,1,'omitnan');
sdD = std(D,0,1,'omitnan');
seD = sdD ./ sqrt(N);
df  = N - 1;
tvals = muD ./ (sdD ./ sqrt(N) + eps);

% ---------------- Thresholds & inference ----------------
switch opt.tail
    case 0,  tthr = tinv(1 - opt.alpha/2, df);
    otherwise, tthr = tinv(1 - opt.alpha, df);
end

switch opt.method
    case "cluster"
        [obsIdx, obsMass] = cluster_masses(tvals, tthr, opt.tail, opt.minClusterPts);

        nullMax = zeros(opt.nPerm,1);
        for p = 1:opt.nPerm
            flips = (rand(N,1)>0.5)*2 - 1; % ±1
            Dp = D .* flips;
            mp = mean(Dp,1); sdp = std(Dp,0,1) + eps;
            tp = mp ./ (sdp ./ sqrt(N));
            [~, massP] = cluster_masses(tp, tthr, opt.tail, opt.minClusterPts);
            nullMax(p) = max([0, massP], [], 'omitnan');
        end
        maxObs = max([0,obsMass], [], 'omitnan');
        p_corr = (sum(nullMax >= maxObs) + 1) / (opt.nPerm + 1);

        Stats = struct('method','cluster','t',tvals,'df',df, ...
                       'tthr',tthr,'clusters',{obsIdx},'massObs',obsMass, ...
                       'maxObs',maxObs,'nullMass',nullMax,'p_corr',p_corr, ...
                       'maskWin',maskWin);

    case "tfce"
        tp = opt.tfceParams;
        if ~isfield(tp,'H'),      tp.H = 2.0; end
        if ~isfield(tp,'E'),      tp.E = 0.5; end
        if ~isfield(tp,'dh'),     tp.dh = 0.1; end
        if ~isfield(tp,'minRun'), tp.minRun = 1; end

        tfceObs = tfce_1d(tvals, opt.tail, tp);

        nullMax = zeros(opt.nPerm,1);
        for p = 1:opt.nPerm
            flips = (rand(N,1)>0.5)*2 - 1;
            Dp = D .* flips;
            mp = mean(Dp,1); sdp = std(Dp,0,1) + eps;
            tpv = mp ./ (sdp ./ sqrt(N));
            sc  = tfce_1d(tpv, opt.tail, tp);
            nullMax(p) = max(sc);
        end
        maxObs = max(tfceObs);
        p_corr = (sum(nullMax >= maxObs) + 1) / (opt.nPerm + 1);
        thr    = quantile(nullMax, 1 - opt.alpha);
        sigMask = tfceObs >= thr;

        Stats = struct('method','tfce','t',tvals,'df',df, ...
                       'tfceScore',tfceObs,'tfceThr',thr,'sigMask',sigMask, ...
                       'nullMax',nullMax,'p_corr',p_corr,'maskWin',maskWin);
    otherwise
        error('Unknown method: %s', opt.method);
end

% ---------------- Plot ----------------
[ax, ~] = ensure_axes(opt.axes); hold(ax,'on'); box(ax,'on');

% (1) Originals (mean±SE)
origBands = gobjects(0); origLines = gobjects(0);
if opt.showOriginal
    % pos
    muPos = mean(TRFpos.W(:,maskWin),1,'omitnan');
    sePos = std(TRFpos.W(:,maskWin),0,1,'omitnan')/sqrt(N);
    y1 = muPos+sePos; y2 = muPos-sePos;
    origBands(end+1) = fill(ax,[tx fliplr(tx)],[y1 fliplr(y2)],opt.posColor, ...
        'FaceAlpha',opt.bandAlpha,'EdgeColor','none','HandleVisibility','off');
    origLines(end+1) = plot(ax,tx,muPos,'-','Color',opt.posColor, ...
        'LineWidth',opt.lineWidth,'DisplayName','pos');

    % neg
    muNeg = mean(TRFneg.W(:,maskWin),1,'omitnan');
    seNeg = std(TRFneg.W(:,maskWin),0,1,'omitnan')/sqrt(N);
    y1 = muNeg+seNeg; y2 = muNeg-seNeg;
    origBands(end+1) = fill(ax,[tx fliplr(tx)],[y1 fliplr(y2)],opt.negColor, ...
        'FaceAlpha',opt.bandAlpha,'EdgeColor','none','HandleVisibility','off');
    origLines(end+1) = plot(ax,tx,muNeg,'-','Color',opt.negColor, ...
        'LineWidth',opt.lineWidth,'DisplayName','neg');
end

% (2) Difference (mean±SE band)
y1 = muD+seD; y2 = muD-seD;
diffBand = fill(ax,[tx fliplr(tx)],[y1 fliplr(y2)],opt.diffColor, ...
    'FaceAlpha',opt.bandAlpha,'EdgeColor','none','HandleVisibility','off');
diffLine = plot(ax,tx,muD,'--','Color',opt.diffColor, ...
    'LineWidth',opt.lineWidth,'DisplayName','diff');

% (3) Zero
if opt.drawZeroLine
    zline = yline(ax,0,':k','HandleVisibility','off');
else
    zline = gobjects(1);
end

% (4) Significant regions
sigPatches = gobjects(0);
[ymin, ymax] = local_yrange(muD,seD, opt.showOriginal, exist('muPos','var'), exist('muNeg','var'), ...
                            exist('sePos','var'), exist('seNeg','var'), ...
                            existvar('muPos'), existvar('sePos'), existvar('muNeg'), existvar('seNeg'));
switch Stats.method
    case 'cluster'
        if ~isempty(Stats.clusters)
            for k = 1:numel(Stats.clusters)
                idx = Stats.clusters{k};
                if isempty(idx), continue; end
                sigPatches(end+1) = patch(ax,[tx(idx) fliplr(tx(idx))], ...
                    [repmat(ymin,1,numel(idx)) repmat(ymax,1,numel(idx))], ...
                    opt.highlightColor,'EdgeColor','none','FaceAlpha',0.25, ...
                    'HandleVisibility','off'); %#ok<AGROW>
            end
        end
    case 'tfce'
        [s,e] = contiguous_runs(Stats.sigMask);
        for k = 1:numel(s)
            idx = s(k):e(k);
            sigPatches(end+1) = patch(ax,[tx(idx) fliplr(tx(idx))], ...
                [repmat(ymin,1,numel(idx)) repmat(ymax,1,numel(idx))], ...
                opt.highlightColor,'EdgeColor','none','FaceAlpha',0.25, ...
                'HandleVisibility','off'); %#ok<AGROW>
        end
end

% (5) Labels/legend
xlabel(ax, opt.xlabel);
ylabel(ax, opt.ylabel);
switch Stats.method
    case 'cluster'
        ttl = sprintf('%s (cluster-corr p_{max}=%.4f)', string(opt.title), Stats.p_corr);
    case 'tfce'
        ttl = sprintf('%s (TFCE: p_{max}=%.4f, thr=%.3g)', string(opt.title), Stats.tfceThr);
end
title(ax, ttl);
if strcmpi(opt.legend,'on')
    legend(ax, [origLines, diffLine], 'Location','best','Box','off','Interpreter','none');
end

H = struct('ax',ax, ...
           'origBands',origBands, 'origLines',origLines, ...
           'diffBand',diffBand, 'diffLine',diffLine, ...
           'zero',zline, 'sig',sigPatches);

end % main

% ===================== helpers =====================

function TRF = select_trf(TRFres, selector)
% Extract one predictor as a unified struct with fields .t (1xT), .W (SxT)
    assert(isstruct(TRFres) && all(isfield(TRFres,{'t','W','predNames'})), ...
        'TRFres must have fields t, W, predNames.');
    t = TRFres.t(:).';
    predNames = string(TRFres.predNames(:)).';
    if isnumeric(selector)
        j = selector;
    else
        j = find(strcmpi(predNames, string(selector)), 1, 'first');
    end
    assert(~isempty(j) && j>=1 && j<=numel(predNames), 'Predictor not found/invalid index.');
    W = squeeze(TRFres.W(:, j, :)); % [S x T]
    % Try to pass mu/se for display if present
    mu = []; se = [];
    if isfield(TRFres,'mu') && ~isempty(TRFres.mu), mu = TRFres.mu(:, j); end
    if isfield(TRFres,'se') && ~isempty(TRFres.se), se = TRFres.se(:, j); end
    TRF = struct('t',t,'W',W,'mu',mu,'se',se,'name',predNames(j));
end

function [ax, fig] = ensure_axes(axIn)
    if isempty(axIn) || ~isvalid(axIn)
        fig = figure('Color','w');
        ax  = axes('Parent',fig);
    else
        ax  = axIn;
        fig = ancestor(ax,'figure');
    end
    hold(ax,'on'); box(ax,'on');
end

function [idxCells, masses] = cluster_masses(tvals, tthr, tail, minPts)
% Supra-threshold clusters and their "mass" (sum of |t| or signed).
    tv = tvals(:).';
    switch tail
        case 0,  supra = abs(tv) >= tthr; weight = abs(tv);
        case 1,  supra = tv >= tthr;      weight = tv;
        case -1, supra = tv <= -tthr;     weight = -tv;
        otherwise, error('tail must be 0, +1, or -1');
    end
    idxCells = {}; masses = [];
    if ~any(supra), return; end
    [s,e] = contiguous_runs(supra);
    for k = 1:numel(s)
        runIdx = s(k):e(k);
        if numel(runIdx) >= minPts
            idxCells{end+1} = runIdx; %#ok<AGROW>
            masses(end+1)   = sum(weight(runIdx), 'omitnan'); %#ok<AGROW>
        end
    end
    if isempty(idxCells), idxCells = {}; masses = []; end
end

function [starts, ends_] = contiguous_runs(mask)
% Start/end indices of true runs in a logical vector (row).
    mask = logical(mask(:).');
    if ~any(mask), starts = []; ends_ = []; return; end
    d = diff([false, mask, false]);
    starts = find(d==1);
    ends_  = find(d==-1) - 1;
end

function score = tfce_1d(tvals, tail, p)
% 1-D TFCE transform of t-values.
    tv = tvals(:).';
    switch tail
        case 0, base = abs(tv);
        case 1, base = max(tv, 0);
        case -1, base = max(-tv, 0);
        otherwise, error('tail must be 0, +1, or -1');
    end
    maxh = max(base);
    if maxh <= 0, score = zeros(size(base)); return; end
    dh = p.dh; H = p.H; E = p.E; minR = p.minRun;
    score = zeros(size(base));
    for h = dh:dh:maxh
        supra = base >= h;
        if ~any(supra), continue; end
        [s,e] = contiguous_runs(supra);
        for k = 1:numel(s)
            runIdx = s(k):e(k);
            if numel(runIdx) < minR, continue; end
            ext = numel(runIdx);
            score(runIdx) = score(runIdx) + (ext^E) * (h^H) * dh;
        end
    end
end

function tf = existvar(v)
% true if variable exists and is non-empty (for local yrange guard)
    tf = exist('v','var') && ~isempty(v); %#ok<EXIST>
end

function [ymin, ymax] = local_yrange(muD,seD, havePosNeg, hasPos, hasNeg, hasSePos, hasSeNeg, muPos, sePos, muNeg, seNeg)
% Safe y-range for significance patches.
    baseMin = min(muD - seD); baseMax = max(muD + seD);
    ymin = baseMin; ymax = baseMax;
    if havePosNeg
        if hasPos && hasSePos && ~isempty(muPos) && ~isempty(sePos)
            ymin = min(ymin, min(muPos - sePos));
            ymax = max(ymax, max(muPos + sePos));
        end
        if hasNeg && hasSeNeg && ~isempty(muNeg) && ~isempty(seNeg)
            ymin = min(ymin, min(muNeg - seNeg));
            ymax = max(ymax, max(muNeg + seNeg));
        end
    end
end
