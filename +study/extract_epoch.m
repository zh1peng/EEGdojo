function Out = extract_study_epoch(varargin)
% EXTRACT_STUDY_EPOCH  Build a single study-level epoch container from many EEGLAB .set files.
%
% Scans for .set files, optionally pre-processes (reference/filter/resample),
% epochs for given markers, applies baseline (if provided), and returns a flat
% container:
%
%   Out.meta   : shared info (fs, times, chanlocs, epoch/baseline window, conditions, preproc, run_tag)
%   Out.<sub>  : fields for each condition only (data matrices, chan × time × trial)
%                e.g., Out.sub_101.x10, Out.sub_101.x20, ...
%
% Name–Value parameters (all optional unless noted):
%   'study_path'     : (char) root dir containing .set files (REQUIRED)
%   'searchstring'   : (char) regex for filenames (default '.*\.set$')
%   'recursive'      : (logical) recurse subfolders (default true)
%   'subject_parser' : (char) regex with named token 'sub' (and optional 'ses') (default '(?<sub>.+)')
%   'markers'        : (cellstr/string) event codes to epoch (REQUIRED)
%   'aliases'        : (2-col cell) {marker, condition}; omitted -> condition == lower(marker)
%   'epoch_window'   : (1x2 double) ms, e.g. [-1000 2000] (default)
%   'baseline'       : (1x2 double or []) ms; [] to skip (default [])
%   'filter_hp'      : (scalar or []) high‑pass Hz; [] skip
%   'filter_lp'      : (scalar or []) low‑pass Hz; [] skip
%   'resample'       : (scalar or []) target Hz; [] keep original
%   'reference'      : (char) 'avg' | 'none' (default 'none')
%   'to_single'      : (logical) cast to single (default false)
%   'save_path'      : (char) optional output .mat path
%   'save_v7_3'      : (logical) save with -v7.3 (default true)
%   'run_tag'        : (char) tag written into Out.meta (default timestamp)
%
% Examples:
%
%   % 1. Minimal example (no preprocessing, no aliases)
%   Out = extract_study_epoch( ...
%       'study_path','/media/NAS/EEGdata/Exp1', ...
%       'searchstring','^sub.*\.set$', ...
%       'markers',{'C101','C102'} );
%
%   % 2. With aliases
%   Out = extract_study_epoch( ...
%       'study_path','/media/NAS/EEGdata/Exp1', ...
%       'markers',{'C101','C102','C103'}, ...
%       'aliases',{'C101','win'; 'C102','loss'; 'C103','neutral'}, ...
%       'epoch_window',[-500 1500], ...
%       'baseline',[-200 0] );
%
%   % 3. With preprocessing (filter + resample + reref)
%   Out = extract_study_epoch( ...
%       'study_path','/media/NAS/EEGdata/Exp2/qced', ...
%       'searchstring','^qced_sub.*\.set$', ...
%       'markers',{'stimA','stimB'}, ...
%       'aliases',{'stimA','target'; 'stimB','nontarget'}, ...
%       'filter_hp',0.1, ...        % high-pass 0.1 Hz
%       'filter_lp',30, ...         % low-pass 30 Hz
%       'resample',250, ...         % downsample to 250 Hz
%       'reference','avg', ...      % average reference
%       'to_single',true, ...
%       'save_path','/media/NAS/EEGdata/Exp2/epochs_all.mat', ...
%       'save_v7_3',true, ...
%       'run_tag','fs250_lp30');
%
%   % 4. With subject/session parsing from filename
%   Out = extract_study_epoch( ...
%       'study_path','/media/NAS/EEGdata/Exp3', ...
%       'searchstring','\.set$', ...
%       'subject_parser','(?<sub>sub\d+)_(?<ses>s\d+)', ...
%       'markers',{'Cue1','Cue2'}, ...
%       'aliases',{'Cue1','left'; 'Cue2','right'});
%
% Requires: EEGLAB on path; your filesearch_regexp on path.

%% ---- Parse Name–Value inputs
p = inputParser; p.FunctionName = 'extract_study_epoch';
addParameter(p,'study_path','',@(x)ischar(x)&&~isempty(x));
addParameter(p,'searchstring','.*\.set$',@ischar);
addParameter(p,'recursive',true,@islogical);
addParameter(p,'subject_parser','(?<sub>.+)',@ischar);

addParameter(p,'markers',{},@(x)iscellstr(x) || isstring(x));
addParameter(p,'aliases',[],@(x)isempty(x) || (iscell(x)&&size(x,2)==2));  % simple 2-col

addParameter(p,'epoch_window',[-1000 2000],@(x)isnumeric(x)&&numel(x)==2); % ms
addParameter(p,'baseline',[],@(x)isnumeric(x) && (isempty(x)||numel(x)==2)); % ms

addParameter(p,'filter_hp',[],@(x)isempty(x)||isscalar(x));
addParameter(p,'filter_lp',[],@(x)isempty(x)||isscalar(x));
addParameter(p,'resample',[],@(x)isempty(x)||isscalar(x));
addParameter(p,'reference','none',@(x)ischar(x)&&~isempty(x));

addParameter(p,'to_single',false,@islogical);
addParameter(p,'save_path','',@ischar);
addParameter(p,'save_v7_3',true,@islogical);
addParameter(p,'run_tag',datestr(now,'yyyy-mm-dd_HHMMSS'),@ischar);

parse(p,varargin{:});
opt = p.Results;

assert(~isempty(opt.study_path)&&isfolder(opt.study_path), 'study_path not found.');
assert(~isempty(opt.markers), 'markers is required.');
markers = cellstr(opt.markers);

% ---- Normalize aliases into two simple parallel lists (no containers.Map)
[alias_mk, alias_cond] = normalizeAliasPairs(opt.aliases);

%% ---- Find .set files (use your filesearch_regexp)
[paths, names] = filesearch_regexp(opt.study_path, opt.searchstring, opt.recursive);
setFiles = fullfile(paths, names);
assert(~isempty(setFiles), 'No .set files matched searchstring in study_path.');

%% ---- Initialize output container
conds = sort(unique(lower(resolveMany(markers, alias_mk, alias_cond))));
Out = struct();
Out.meta = struct( ...
    'fs',              [], ...
    'times',           [], ...
    'chanlocs',        [], ...
    'epoch_window_ms', opt.epoch_window, ...
    'baseline_ms',     opt.baseline, ...
    'preproc',         struct('filter_hp',opt.filter_hp,'filter_lp',opt.filter_lp, ...
                              'resample',opt.resample,'reference',opt.reference), ...
    'created_at',      datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
    'run_tag',         opt.run_tag, ...
    'version',         1);
Out.meta.conditions = conds;

%% ---- Iterate files (per subject)
firstMetaLocked = false;
warns = {};

% summary collectors
summary_sub  = {};
summary_cond = {};
summary_n    = [];
summary_file = {};

for i = 1:numel(setFiles)
    fpath = setFiles{i};
    [subID, ~] = parseSubSesFromFile(fpath, opt.subject_parser);
    subKey = makeFieldKey(subID);

    if ~isfield(Out, subKey)
        Out.(subKey) = struct();   % no subject_id/session fields
    end

    % Load
    EEG = pop_loadset(fpath);
    EEG = eeg_checkset(EEG);

    % ---- Reference (inline)
    if ~isempty(opt.reference) && ~strcmpi(opt.reference,'none')
        switch lower(opt.reference)
            case 'avg'
                EEG = pop_reref(EEG, []);
            otherwise
                warning('Unknown reference mode "%s" (using none).', opt.reference);
        end
        EEG = eeg_checkset(EEG);
    end

    % ---- Filter (inline)
    if ( ~isempty(opt.filter_hp) && opt.filter_hp>0 ) || ( ~isempty(opt.filter_lp) && opt.filter_lp>0 )
        args = {};
        if ~isempty(opt.filter_hp) && opt.filter_hp>0, args = [args, {'locutoff', opt.filter_hp}]; end
        if ~isempty(opt.filter_lp) && opt.filter_lp>0, args = [args, {'hicutoff', opt.filter_lp}]; end
        EEG = pop_eegfiltnew(EEG, args{:});
        EEG = eeg_checkset(EEG);
    end

    % ---- Resample (inline)
    if ~isempty(opt.resample) && opt.resample>0 && EEG.srate~=opt.resample
        EEG = pop_resample(EEG, opt.resample);
        EEG = eeg_checkset(EEG);
    end

    % Lock meta from the first processed file
    if ~firstMetaLocked
        Out.meta.fs       = EEG.srate;
        Out.meta.chanlocs = EEG.chanlocs;
        firstMetaLocked   = true;
    else
        if EEG.srate ~= Out.meta.fs
            warns{end+1} = sprintf('[WARN] %s: srate %g != meta.fs %g', subID, EEG.srate, Out.meta.fs);
        end
        if numel(EEG.chanlocs) ~= numel(Out.meta.chanlocs)
            warns{end+1} = sprintf('[WARN] %s: channel count differs from first file', subID);
        end
    end

    % Per-marker epoch loop
    for m = 1:numel(markers)
        mk       = markers{m};
        condName = resolveAlias(mk, alias_mk, alias_cond);
        condKey  = makeFieldKey(condName);

        try
            EEGep = pop_epoch(EEG, {mk}, opt.epoch_window/1000);  % ms -> s
            EEGep = eeg_checkset(EEGep);

            if ~isempty(opt.baseline)
                EEGep = pop_rmbase(EEGep, opt.baseline, []);      % ms
            end

            if isempty(Out.meta.times) && isfield(EEGep,'times')
                Out.meta.times = EEGep.times;
            end

            data = EEGep.data;                                    % chan × time × trial
            if opt.to_single, data = single(data); end

            % ---- FLAT storage: only data
            Out.(subKey).(condKey) = data;


            % ---- add to summary
            summary_sub{end+1,1}  = subID; 
            summary_cond{end+1,1} = condKey;
            summary_n(end+1,1)    = size(data,3);
            summary_file{end+1,1} = fpath;

        catch ME
            warns{end+1} = sprintf('[WARN] %s / %s: %s', subID, condName, ME.message);
            Out.(subKey).(condKey) = [];
            % add zero-row to summary
            summary_sub{end+1,1}  = subID;
            summary_cond{end+1,1} = condKey;
            summary_n(end+1,1)    = 0;
            summary_file{end+1,1} = fpath;
        end
    end
end


% ==== build WIDE summary table (one row per subject; cols = conditions) ====
% Tall table first
T = table(summary_sub, summary_cond, summary_n, summary_file, ...
    'VariableNames', {'sub','condition','n','file'});

% Pivot trial counts to wide by condition
Tw = unstack(T(:, {'sub','condition','n'}), 'n', 'condition', 'GroupingVariables', 'sub');

% Ensure all expected condition columns exist (even if some subjects lack them)
allConds = conds(:)';  % from earlier
for k = 1:numel(allConds)
    vn = allConds{k};
    if ~ismember(vn, Tw.Properties.VariableNames)
        Tw.(vn) = zeros(height(Tw),1);
    end
end

% Order columns: sub, <conds...>
Tw = Tw(:, ['sub', allConds]);

% One file path per subject (first seen) — robust across MATLAB versions
[uniqSubs, ia] = unique(T.sub, 'stable');
firstFiles     = T.file(ia);
F = table(uniqSubs, firstFiles, 'VariableNames', {'sub','file'});

% Join files onto wide table
Tw = outerjoin(Tw, F, 'Keys', 'sub', 'MergeKeys', true);

% Replace missing trial counts with 0 (defensive)
Tw{:, allConds} = fillmissing(Tw{:, allConds}, 'constant', 0);

% Attach to output
Out.meta.trialN = Tw;

%% ---- Save if requested
if ~isempty(opt.save_path)
    if opt.save_v7_3
        save(opt.save_path, 'Out', '-v7.3');
    else
        save(opt.save_path, 'Out');
    end
end

% Optional: print warnings
if ~isempty(warns)
    fprintf('=== Extraction warnings (%d) ===\n', numel(warns));
    fprintf('%s\n', strjoin(warns,newline));
end

end  % main function


%% ----------------- Local helpers (kept minimal & readable) -----------------

function [alias_mk, alias_cond] = normalizeAliasPairs(aliases)
% Accepts [] or a 2-col cell {marker, condition}; returns two cell arrays.
    if isempty(aliases)
        alias_mk = {}; alias_cond = {};
        return;
    end
    if iscell(aliases) && size(aliases,2)==2
        alias_mk   = cellstr(aliases(:,1));
        alias_cond = cellstr(aliases(:,2));
        return;
    end
    error('aliases must be a 2-col cell {marker, condition}.');
end

function out = resolveMany(markers, alias_mk, alias_cond)
    out = cell(size(markers));
    for i=1:numel(markers)
        out{i} = resolveAlias(markers{i}, alias_mk, alias_cond);
    end
end

function name = resolveAlias(marker, alias_mk, alias_cond)
% If marker exists in alias list (case-insensitive), use the mapped condition; else lower(marker).
    if ~isempty(alias_mk)
        idx = find(strcmpi(marker, alias_mk), 1, 'first');
        if ~isempty(idx)
            name = alias_cond{idx};
            name = lower(name);
            return;
        end
    end
    name = lower(marker);
end

function [subID, sesID] = parseSubSesFromFile(fpath, parserPattern)
    [~,fname,~] = fileparts(fpath);
    m = regexp(fname, parserPattern, 'names');
    assert(~isempty(m) && isfield(m,'sub'), 'Failed to parse subject from: %s', fname);
    subID = m.sub;
    if isfield(m,'ses') && ~isempty(m.ses), sesID = m.ses; else, sesID = ''; end
end

function key = makeFieldKey(txt)
    key = lower(regexprep(txt,'[^A-Za-z0-9]','_'));
    if ~isempty(key) && isstrprop(key(1),'digit'), key = ['x' key]; end
end
