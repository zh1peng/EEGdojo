function [summary_table, analysis_obj] = explore(analysis_obj, varargin)
% EXPLORE  Batch-run contrasts, stats, and plots for ERPanalysis/TFanalysis.
%
% [summary_table, analysis_obj] = study.explore(analysis_obj, Name,Value,...)
%
% Works with:
%   - study.ERPanalysis  (ERP: compute_erps → compute_ga → compute_contrast)
%   - study.TFanalysis   (TF : compute_ga_tfd → define_contrast)
%
% Name-Value options
% -------------------------------------------------------------------------
% 'output_dir'   (char, required)  Where to save figures and CSV summary.
% 'stat_method'  (function_handle)  Default auto-selected:
%                   ERPanalysis -> @compute_stats
%                   TFanalysis  -> @compute_tfce_stats (override if desired)
% 'stat_params'  (struct)        Name/value pairs for the stats function.
%                                e.g., struct('alpha',0.05,'mcc','fdr','roi','Cz_ROI')
%
% 'plot_method'  (function_handle)  Default auto-selected:
%                   ERPanalysis -> @plot_contrast_erp
%                   TFanalysis  -> @plot_contrast_tfr
% 'plot_target'  (char, required)   Channel or ROI name (e.g., 'Cz' or 'Midline').
% 'plot_params'  (cell)             Extra plotting NV pairs, e.g., {'show_sig',true}
%
% 'include_within'  (logical)   Default true.  All condition pairs within each group.
% 'include_between' (logical)   Default true.  All group pairs within each condition.
% 'save_format'     (char)      'png'|'pdf'|'tif' (default 'png')
% 'dpi'             (double)    Raster resolution (default 150)
% 'only_significant_plots' (logical) Default true.  Save figure only if significant.
%
% Returns
% -------------------------------------------------------------------------
% summary_table  Table with columns:
%   ContrastName, PositiveTerm, NegativeTerm, IsPaired, Npos, Nneg,
%   StatType, SignificantCount, FigureFile
%
% Notes
% -------------------------------------------------------------------------
% • For ERPanalysis, plotting expects 'target_name' as a name/value.
%   For TFanalysis, plotting expects the target as a positional argument.
%   This function tries positional first, then falls back to name/value.
% • The function mutates analysis_obj (adds contrasts and stats).
%
% Example
% -------------------------------------------------------------------------
% erp = study.ERPanalysis(ds).define_group('All', ds.get_subjects()) ...
%     .select_conditions({'A','B'}).compute_erps().compute_ga();
% stat_p = struct('alpha',0.05,'mcc','fdr','roi','Cz_ROI');
% plot_p = {'show_sig',true,'show_diff',true};
% [T, erp] = study.explore(erp, 'output_dir','./res', ...
%     'stat_params',stat_p, 'plot_params',plot_p, 'plot_target','Cz');

% ------------ Parse inputs ------------
p = inputParser;
addRequired(p, 'analysis_obj');
addParameter(p, 'output_dir', '', @ischar);
addParameter(p, 'stat_method', [], @(f) isempty(f) || isa(f,'function_handle'));
addParameter(p, 'stat_params', struct(), @isstruct);
addParameter(p, 'plot_method', [], @(f) isempty(f) || isa(f,'function_handle'));
addParameter(p, 'plot_target', '', @ischar);
addParameter(p, 'plot_params', {}, @iscell);
addParameter(p, 'include_within', true, @islogical);
addParameter(p, 'include_between', true, @islogical);
addParameter(p, 'save_format', 'png', @(s)ischar(s)&&ismember(lower(s),{'png','pdf','tif','tiff'}));
addParameter(p, 'dpi', 150, @isnumeric);
addParameter(p, 'only_significant_plots', true, @islogical);
parse(p, analysis_obj, varargin{:});
opt = p.Results;

if isempty(opt.output_dir), error('Output directory is required.'); end
if ~isfolder(opt.output_dir), mkdir(opt.output_dir); end
if isempty(opt.plot_target), error('''plot_target'' (channel/ROI) is required.'); end

% ------------ Auto-detect class & defaults ------------
cls = class(analysis_obj);
isERP = contains(cls, 'ERPanalysis');
isTF  = contains(cls, 'TFanalysis');
if ~(isERP || isTF)
    error('analysis_obj must be ERPanalysis or TFanalysis.');
end

% Default methods if not provided
if isempty(opt.stat_method)
    opt.stat_method = iff(isERP, @compute_stats, @compute_tfce_stats);
end
if isempty(opt.plot_method)
    opt.plot_method = iff(isERP, @plot_contrast_erp, @plot_contrast_tfr);
end

% ------------ Ensure GA exists (and prerequisites) ------------
if isERP
    if ~isfield(analysis_obj.Results,'ERPs') || isempty(fieldnames(analysis_obj.Results.ERPs))
        warning('Subject ERPs missing. Calling compute_erps()...');
        analysis_obj = analysis_obj.compute_erps();
    end
    if ~isfield(analysis_obj.Results,'GA') || isempty(fieldnames(analysis_obj.Results.GA))
        warning('Grand Averages missing. Calling compute_ga()...');
        analysis_obj = analysis_obj.compute_ga();
    end
else % TF
    if ~isfield(analysis_obj.Results,'GA_TFD') || isempty(fieldnames(analysis_obj.Results.GA_TFD))
        warning('Grand-average TFRs missing. Calling compute_ga_tfd()...');
        analysis_obj = analysis_obj.compute_ga_tfd();
    end
end

% ------------ Build contrast list ------------
groups = fieldnames(analysis_obj.Selection.Groups);
conditions = analysis_obj.Selection.Conditions;
contrast_list = {};

% Within-group: all pairs of conditions
if opt.include_within && numel(conditions) >= 2
    pairs = nchoosek(1:numel(conditions), 2);
    for g = 1:numel(groups)
        for k = 1:size(pairs,1)
            c1 = conditions{pairs(k,1)}; c2 = conditions{pairs(k,2)};
            cx.name      = sprintf('%s__%s_vs_%s', groups{g}, c1, c2);
            cx.positive  = {groups{g}, c1};
            cx.negative  = {groups{g}, c2};
            cx.is_paired = true;  % same subjects by definition of group
            contrast_list{end+1} = cx; %#ok<AGROW>
        end
    end
end

% Between-group: all pairs of groups per condition
if opt.include_between && numel(groups) >= 2
    gpairs = nchoosek(1:numel(groups), 2);
    for c = 1:numel(conditions)
        for k = 1:size(gpairs,1)
            g1 = groups{gpairs(k,1)}; g2 = groups{gpairs(k,2)};
            cx.name      = sprintf('%s_vs_%s__%s', g1, g2, conditions{c});
            cx.positive  = {g1, conditions{c}};
            cx.negative  = {g2, conditions{c}};
            % paired only if subject lists are identical
            s1 = analysis_obj.Selection.Groups.(g1);
            s2 = analysis_obj.Selection.Groups.(g2);
            cx.is_paired = isequal(s1, s2);
            contrast_list{end+1} = cx; %#ok<AGROW>
        end
    end
end

fprintf('Generated %d contrasts to explore.\n', numel(contrast_list));
if isempty(contrast_list)
    summary_table = table();
    return;
end

% ------------ Helper: struct -> name/value cell ------------
stat_args = struct2nvpairs(opt.stat_params);

% ------------ Figure save helper ------------
function save_figure_safe(fig, path, fmt, dpi)
    [~,~,ext] = fileparts(path);
    if isempty(ext), path = [path '.' fmt]; end
    try
        if exist('exportgraphics','file') && any(strcmpi(fmt,{'png','pdf','tif','tiff'}))
            exportgraphics(fig, path, 'Resolution', dpi);
        else
            switch lower(fmt)
                case 'png',  print(fig, path, '-dpng',  ['-r' num2str(dpi)]);
                case {'tif','tiff'}, print(fig, path, '-dtiff',['-r' num2str(dpi)]);
                case 'pdf',  print(fig, path, '-dpdf');
            end
        end
    catch ME
        warning('export failed (%s); falling back to saveas: %s', fmt, ME.message);
        saveas(fig, path);
    end
end

% ------------ Contrast runner ------------
summary_rows = {};
for i = 1:numel(contrast_list)
    cx = contrast_list{i};
    fprintf('---\nRunning contrast: %s\n', cx.name);

    % 1) Define/compute contrast on the analysis object
    if any(strcmp(methods(analysis_obj), 'compute_contrast'))
        analysis_obj = analysis_obj.compute_contrast(cx.name, cx.positive, cx.negative);
    elseif any(strcmp(methods(analysis_obj), 'define_contrast'))
        analysis_obj = analysis_obj.define_contrast(cx.name, cx.positive, cx.negative);
    else
        error('Neither compute_contrast nor define_contrast found on analysis_obj.');
    end

    % 2) Stats
    try
        analysis_obj = feval(opt.stat_method, analysis_obj, 'contrast', cx.name, stat_args{:});
        stat_type = func2str(opt.stat_method);
    catch ME
        warning('Stats failed for %s: %s', cx.name, ME.message);
        stat_type = [func2str(opt.stat_method) ' (FAILED)'];
    end

    % 3) Significance check
    sig_count = 0; is_sig = false;
    figsav = 'N/A';
    res = analysis_obj.Results.Contrasts.(cx.name);

    if isfield(res, 'Stats_TFCE') && isstruct(res.Stats_TFCE) && isfield(res.Stats_TFCE,'cluster_p')
        cp = res.Stats_TFCE.cluster_p;
        a  = res.Stats_TFCE.alpha;
        if isempty(a), a = 0.05; end
        sig_idx = find(cp < a);
        sig_count = numel(sig_idx);
        is_sig = sig_count > 0;
    elseif isfield(res,'Stats') && isstruct(res.Stats) && isfield(res.Stats,'h')
        sig_count = sum(res.Stats.h(:));
        is_sig = sig_count > 0;
    end

    % 4) Plot (only if significant, unless overridden)
    if (~opt.only_significant_plots) || is_sig
        try
            % Try positional target (TFanalysis style)
            fig = figure('Visible','off'); close(fig); % ensure gcf resets later
            feval(opt.plot_method, analysis_obj, cx.name, opt.plot_target, opt.plot_params{:});
        catch
            % Fallback: name/value target (ERPanalysis style)
            feval(opt.plot_method, analysis_obj, cx.name, 'target_name', opt.plot_target, opt.plot_params{:});
        end
        fig = gcf; set(fig,'Color','w');
        safe_name = regexprep(cx.name,'[^a-zA-Z0-9_\-]','_');
        fpath = fullfile(opt.output_dir, sprintf('%s.%s', safe_name, lower(opt.save_format)));
        save_figure_safe(fig, fpath, lower(opt.save_format), opt.dpi);
        close(fig);
        figsav = fpath;
    end

    % 5) Npos/Nneg (best-effort)
    Npos = NaN; Nneg = NaN;
    if isERP
        try
            Npos = analysis_obj.Results.GA.(cx.positive{1}).(cx.positive{2}).n;
            Nneg = analysis_obj.Results.GA.(cx.negative{1}).(cx.negative{2}).n;
        catch, end
    else
        if isfield(res,'n_pos'), Npos = res.n_pos; end
        if isfield(res,'n_neg'), Nneg = res.n_neg; end
    end

    summary_rows(end+1,:) = {cx.name, sprintf('%s:%s',cx.positive{1},cx.positive{2}), ...
        sprintf('%s:%s',cx.negative{1},cx.negative{2}), cx.is_paired, Npos, Nneg, ...
        stat_type, sig_count, figsav}; %#ok<AGROW>
end

summary_table = cell2table(summary_rows, 'VariableNames', ...
    {'ContrastName','PositiveTerm','NegativeTerm','IsPaired','Npos','Nneg','StatType','SignificantCount','FigureFile'});

% Save CSV
csv_path = fullfile(opt.output_dir, 'exploration_summary.csv');
writetable(summary_table, csv_path);
fprintf('---\nExploration complete. Summary saved to %s\n', opt.output_dir);

end % function explore

% ---------------- local helpers ----------------
function out = iff(cond,a,b), if cond, out=a; else, out=b; end
end

function nv = struct2nvpairs(s)
% Flatten a scalar struct into name/value cell array.
if isempty(s)
    nv = {};
    return;
end
fn = fieldnames(s);
vals = struct2cell(s);
nv = reshape([fn.'; vals.'], 1, []); % interleave name/value
end
