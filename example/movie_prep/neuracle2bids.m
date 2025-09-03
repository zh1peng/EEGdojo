function neuracle2bids(rootSourcedata, params)
% Convert Neuracle sourcedata to BIDS-style EEGLAB .set + sidecars (RAW tree).
%
% Required:
%   rootSourcedata            e.g., 'V:\BEAT_BIDS\sourcedata'
%   params.ChanLocation       path to .elc/.ced (e.g., standard_1005.elc)
%
% Optional:
%   params.OutRoot            BIDS RAW root (default = parent(rootSourcedata))
%   params.WriteSidecars      true/false (default = true)
%
% Notes:
% - Finds folders containing 'data.bdf' (and optional 'evt.bdf')
% - Derives BIDS names from path: ...\sourcedata\sub-101\1_rest -> sub-101_task-rest_run-01
% - If impedance JSON is found in the task folder, values are written to _electrodes.tsv (kΩ)
% - If no impedance, filenames are unchanged; sidecars record absence

if nargin < 2 || ~isfield(params,'ChanLocation')
    error('params.ChanLocation is required (path to standard_1005.elc or .ced)');
end
if ~isfield(params,'OutRoot') || isempty(params.OutRoot)
    % RAW BIDS root (sibling of sourcedata)
    params.OutRoot = fileparts(rootSourcedata);
end
if ~exist(params.OutRoot,'dir'); mkdir(params.OutRoot); end
if ~isfield(params,'WriteSidecars'); params.WriteSidecars = true; end

summaryRows = [];

recFolders = find_recording_folders(rootSourcedata);
fprintf('Found %d recording folder(s).\n', numel(recFolders));
for i = 1:numel(recFolders)
    inputPath = recFolders{i};
    try
        [bidsStem, outDir, sub, taskLabel, runNum] = derive_bids_paths(inputPath, params.OutRoot);

        % --- Import & capture importer setname BEFORE changing it
        [EEG, hadEvents] = import_neuracle(inputPath);
        importSetname = EEG.setname;

        % Assign channel locations
        EEG = pop_chanedit(EEG, 'lookup', params.ChanLocation);
        EEG = eeg_checkset(EEG);
        EEG.urchanlocs = EEG.chanlocs;

        % Read impedance map if present
        [impMap, hadImpedance, impOnlineFlag] = read_impedance_map(inputPath);

        % Set a clean setname (BIDS stem)
        EEG.setname = bidsStem;

        % Ensure output folder exists (.../sub-xxx/eeg)
        if ~exist(outDir,'dir'); mkdir(outDir); end

        % Save .set (RAW BIDS tree)
        setFile = [bidsStem '_eeg.set'];
        EEG = pop_saveset(EEG, 'filename', setFile, 'filepath', outDir, 'check','on');
        outSetPath = fullfile(outDir, setFile);
        fprintf('Saved: %s\n', outSetPath);

        % Sidecars
        if params.WriteSidecars
            % write_events_tsv(EEG, fullfile(outDir, [bidsStem '_events.tsv']), hadEvents);
            write_channels_tsv(EEG, fullfile(outDir, [bidsStem '_channels.tsv']));
            write_electrodes_and_coords(EEG, ...
                fullfile(outDir, [bidsStem '_electrodes.tsv']), ...
                fullfile(outDir, [bidsStem '_coordsystem.json']), impMap, sub);

            % NEW: impedance meta (present/absent + source hint)
            write_impedance_json(fullfile(outDir, [bidsStem '_impedance.json']), ...
                                 hadImpedance, impOnlineFlag);
        end

        % Sanity check: setname vs parsed subject
        if ~contains(lower(importSetname), lower(sub))
            warning('Setname/subject mismatch: import setname "%s" vs parsed subject "%s" (%s)', ...
                    importSetname, sub, inputPath);
        end

        % Append summary row
        durMin = EEG.pnts / EEG.srate / 60;
        summaryRows = [summaryRows; struct( ...
            'subject', string(sub), ...
            'task',    string(taskLabel), ...
            'run',     runNum, ...
            'bids_stem', string(bidsStem), ...
            'import_setname', string(importSetname), ...
            'nbchan', EEG.nbchan, ...
            'srate',  EEG.srate, ...
            'duration_min', durMin, ...
            'n_events', numel(EEG.event), ...
            'impedance_present', logical(hadImpedance), ...
            'input_path', string(inputPath), ...
            'output_set', string(outSetPath) ...
        )];

    catch ME
        warning('Failed on %s: %s', inputPath, ME.message);
    end
end

% Summary TSV (under RAW root)
if ~isempty(summaryRows)
    T = struct2table(summaryRows);
    summaryFile = fullfile(params.OutRoot, 'conversion_summary.tsv');
    writetable(T, summaryFile, 'FileType','text','Delimiter','\t');
    fprintf('Wrote summary: %s\n', summaryFile);

    % Duplicate check: subject+task+run
    [grp,~,idx] = unique(T(:,{'subject','task','run'}), 'rows');
    counts = accumarray(idx, 1);
    dups = grp(counts > 1, :);
    if ~isempty(dups)
        fprintf('WARNING: Duplicate subject-task-run combinations found:\n');
        disp(dups);
    end
end
end

% -------- Helpers --------

function recFolders = find_recording_folders(rootSourcedata)
D = dir(fullfile(rootSourcedata, 'sub-*'));
recFolders = {};
for s = 1:numel(D)
    if ~D(s).isdir, continue; end
    subdir = fullfile(D(s).folder, D(s).name);
    T = dir(fullfile(subdir, '*'));
    for t = 1:numel(T)
        if T(t).isdir
            cand = fullfile(T(t).folder, T(t).name);
            if exist(fullfile(cand, 'data.bdf'),'file')
                recFolders{end+1} = cand; %#ok<AGROW>
            end
        end
    end
end
end
function [bidsStem, outDir, sub, taskLabel, runNum] = derive_bids_paths(inputPath, outRoot)
parts = strsplit(inputPath, filesep);
subIdx = find(startsWith(parts, 'sub-'), 1, 'last');
if isempty(subIdx); error('No sub-* in path: %s', inputPath); end
sub = parts{subIdx};
leaf = parts{end};

% Parse leading run number and the rest as raw task label (allow underscores etc.)
% Matches "5_Movie_P", "6-Movie-DM", "7 movie dm", "12rest", etc.
tok = regexp(leaf, '^(\d+)[\s\-_]*(.*)$', 'tokens', 'once');
if ~isempty(tok)
    runNum = str2double(tok{1});
    rawTask = strtrim(tok{2});
    if isempty(rawTask)
        rawTask = 'rest';
    end
else
    % Fallback: no leading number; treat whole leaf as task, run=1
    runNum = 1;
    rawTask = leaf;
end

% BIDS task label must be alphanumeric only; strip everything else and lowercase
taskLabel = lower(regexprep(rawTask, '[^A-Za-z0-9]+', ''));

% Final safety fallback
if isempty(taskLabel)
    taskLabel = 'rest';
end

% here we set run number as 1 per task
runNum=1;

runStr  = sprintf('run-%02d', runNum);
taskStr = sprintf('task-%s', taskLabel);
bidsStem = sprintf('%s_%s_%s', sub, taskStr, runStr);

% RAW BIDS output dir: <BIDS root>/sub-xxx/eeg
outDir = fullfile(outRoot, sub, 'eeg');
end


function [EEG, hadEvents] = import_neuracle(inputPath)
dataFile = fullfile(inputPath, 'data.bdf');
evtFile  = fullfile(inputPath, 'evt.bdf');
if ~exist(dataFile,'file'); error('data.bdf not found in %s', inputPath); end
filenames = {'data.bdf'};
if exist(evtFile,'file'); filenames = {'data.bdf','evt.bdf'}; end
EEG = pop_importNeuracle(filenames, inputPath, ...
    'channels', [], 'blockrange', [], ...
    'importevent', 'on', 'memorymapped', 'off');
EEG = eeg_checkset(EEG);
EEG = eeg_checkset(EEG, 'eventconsistency');
hadEvents = ~isempty(EEG.event);
end

function [impMap, had, impOnlineFlag] = read_impedance_map(taskFolder)
% Also tries to read "ImpedanceOnline" if present in the JSON.
impMap = containers.Map('KeyType','char','ValueType','double');
had = false; impOnlineFlag = [];
J = dir(fullfile(taskFolder, '*.json'));
for k = 1:numel(J)
    try
        txt = fileread(fullfile(taskFolder, J(k).name));
        S = jsondecode(txt);
        if isfield(S,'ImpedanceOnline'); impOnlineFlag = logical(S.ImpedanceOnline); end
        if isfield(S,'ImpedanceData') && isfield(S.ImpedanceData,'Item1') && isfield(S.ImpedanceData,'Item2')
            labs = string(S.ImpedanceData.Item1);
            vals = double(S.ImpedanceData.Item2);
            if numel(labs)==numel(vals)
                for i=1:numel(labs)
                    key = lower(char(labs(i)));
                    val = vals(i);
                    if isfinite(val) && val>=0 && val<1e6
                        impMap(key) = val; % assumed kΩ
                    end
                end
                had = ~isempty(impMap.keys);
                if had, return; end
            end
        end
    catch
        % ignore and continue
    end
end
end

% function write_events_tsv(EEG, outFile, hadEvents)
% if ~hadEvents, return; end
% onset    = arrayfun(@(e) (e.latency-1)/EEG.srate, EEG.event)';
% hasDur   = arrayfun(@(e) isfield(e,'duration') && ~isempty(e.duration), EEG.event);
% duration = nan(numel(EEG.event),1);
% duration(hasDur) = [EEG.event(hasDur).duration]'/EEG.srate;
% trial_type = string({EEG.event.type})';
% T = table(onset, duration, trial_type);
% writetable(T, outFile, 'FileType','text','Delimiter','\t');
% end

function write_channels_tsv(EEG, outFile)
names = string({EEG.chanlocs.labels})';
types = strings(EEG.nbchan,1);
for i=1:EEG.nbchan
    if isfield(EEG.chanlocs(i),'type') && ~isempty(EEG.chanlocs(i).type)
        types(i) = string(EEG.chanlocs(i).type);
    else
        types(i) = "EEG";
    end
end
status = repmat("good", EEG.nbchan,1);
Tch = table(names, types, status);
writetable(Tch, outFile, 'FileType','text','Delimiter','\t');
end

function write_electrodes_and_coords(EEG, electrodesFile, coordsFile, impMap, sub)
names = string({EEG.chanlocs.labels})';
hasXYZ = isfield(EEG.chanlocs,'X') && ~isempty(EEG.chanlocs(1).X);

% electrodes.tsv (with optional impedance in kΩ, per BIDS)
if hasXYZ
    X = arrayfun(@(c) c.X, EEG.chanlocs)'; Y = arrayfun(@(c) c.Y, EEG.chanlocs)';
    Z = arrayfun(@(c) c.Z, EEG.chanlocs)';
else
    X = nan(EEG.nbchan,1); Y = X; Z = X;
end

impedance = nan(EEG.nbchan,1);
if ~isempty(impMap)
    lblLower = lower(string({EEG.chanlocs.labels}));
    for i=1:EEG.nbchan
        key = char(lblLower(i));
        if isKey(impMap, key)
            impedance(i) = impMap(key);  % kΩ (per your JSON)
        end
    end
end

Tel = table(names, X, Y, Z, impedance);
writetable(Tel, electrodesFile, 'FileType','text','Delimiter','\t');

% coordsystem.json (minimal)
cs = struct( ...
    'EEGCoordinateSystem', 'EEGLAB', ...
    'EEGCoordinateUnits',  'mm', ...
    'EEGCoordinateSystemDescription', ...
      'EEGLAB DIPFIT standard_1005.elc (template; BEM/MNI-based)', ...
    'InstitutionSubject', sub);
fid = fopen(coordsFile, 'w');
fwrite(fid, jsonencode(cs, 'PrettyPrint', true));
fclose(fid);
end

function write_impedance_json(outFile, hadImpedance, impOnlineFlag)
S = struct('ImpedancePresent', logical(hadImpedance));
if ~isempty(impOnlineFlag)
    S.ImpedanceOnline = logical(impOnlineFlag);
end
fid = fopen(outFile, 'w');
fwrite(fid, jsonencode(S, 'PrettyPrint', true));
fclose(fid);
end
