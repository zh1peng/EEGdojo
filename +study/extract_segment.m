function Out = extract_segment(varargin)
% EXTRACT_SEGMENT Build a study-level container from cropped segments of many EEGLAB .set files.
% 
% This function scans for .set files, applies `prep.crop_by_markers` to extract a
% continuous segment based on specified markers, and returns a flat container.
% 
%   Out.meta   : shared info (fs, srate, chanlocs, markers, padding, etc.)
%   Out.<sub>  : fields for each subject (e.g., Out.sub_101) containing the
%                cropped data matrix (chan × time).
% 
% Name–Value parameters (all optional unless noted):
%   'study_path'     : (char) root dir containing .set files (REQUIRED)
%   'searchstring'   : (char) regex for filenames (default '.*\.set$')
%   'recursive'      : (logical) recurse subfolders (default true)
%   'subject_parser' : (char) regex with named token 'sub' (default '(?<sub>.+)')
%   'StartMarker'    : (char | string) The marker for the start of the segment (REQUIRED).
%   'EndMarker'      : (char | string) The marker for the end of the segment (REQUIRED).
%   'PadSec'         : (numeric) Seconds to pad on both sides (default 0).
%   'to_single'      : (logical) cast to single (default false)
%   'save_path'      : (char) optional output .mat path
%   'save_v7_3'      : (logical) save with -v7.3 (default true)
%   'run_tag'        : (char) tag written into Out.meta (default timestamp)
% 
% Example:
%   Out = study.extract_segment( ...
%       'study_path', 'Z:\path\to\your\data', ...
%       'StartMarker', 'movie_start',
%       'EndMarker', 'movie_end',
%       'PadSec', 5);
% 
% Requires: EEGLAB on path; filesearch_regexp and prep.crop_by_markers on path.

%% ---- Parse Name–Value inputs
p = inputParser;
p.FunctionName = 'extract_segment';
addParameter(p, 'study_path', '', @(x) ischar(x) && ~isempty(x));
addParameter(p, 'searchstring', '.*\.set$', @ischar);
addParameter(p, 'recursive', true, @islogical);
addParameter(p, 'subject_parser', '(?<sub>.+)', @ischar);

addParameter(p, 'StartMarker', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'EndMarker', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'PadSec', 0, @isnumeric);

addParameter(p, 'to_single', false, @islogical);
addParameter(p, 'save_path', '', @ischar);
addParameter(p, 'save_v7_3', true, @islogical);
addParameter(p, 'run_tag', datestr(now, 'yyyy-mm-dd_HHMMSS'), @ischar);

parse(p, varargin{:});
opt = p.Results;

assert(~isempty(opt.study_path) && isfolder(opt.study_path), 'study_path not found.');
assert(~isempty(opt.StartMarker) && ~isempty(opt.EndMarker), 'StartMarker and EndMarker are required.');

%% ---- Find .set files
[paths, names] = filesearch_regexp(opt.study_path, opt.searchstring, opt.recursive);
setFiles = fullfile(paths, names);
assert(~isempty(setFiles), 'No .set files matched searchstring in study_path.');

%% ---- Initialize output container
Out = struct();
Out.meta = struct( ...
    'fs',              [], ...
    'chanlocs',        [], ...
    'start_marker',    opt.StartMarker, ...
    'end_marker',      opt.EndMarker, ...
    'padding_sec',     opt.PadSec, ...
    'created_at',      datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
    'run_tag',         opt.run_tag, ...
    'version',         1);

%% ---- Iterate files (per subject)
firstMetaLocked = false;
warns = {};

for i = 1:numel(setFiles)
    fpath = setFiles{i};
    [subID, ~] = parseSubSesFromFile(fpath, opt.subject_parser);
    subKey = makeFieldKey(subID);

    try
        % Load
        EEG = pop_loadset(fpath);
        EEG = eeg_checkset(EEG);

        % Crop the segment
        [EEG_cropped, crop_info] = prep.crop_by_markers(EEG, ...
            'StartMarker', opt.StartMarker, ...
            'EndMarker', opt.EndMarker, ...
            'PadSec', opt.PadSec);

        % Lock meta from the first processed file
        if ~firstMetaLocked
            Out.meta.fs = EEG_cropped.srate;
            Out.meta.chanlocs = EEG_cropped.chanlocs;
            firstMetaLocked = true;
        else
            if EEG_cropped.srate ~= Out.meta.fs
                warns{end+1} = sprintf('[WARN] %s: srate %g ~= meta.fs %g', subID, EEG_cropped.srate, Out.meta.fs);
            end
            if numel(EEG_cropped.chanlocs) ~= numel(Out.meta.chanlocs)
                warns{end+1} = sprintf('[WARN] %s: channel count differs from first file', subID);
            end
        end

        data = EEG_cropped.data;
        if opt.to_single, data = single(data); end

        % Store data
        Out.(subKey).data = data;
        Out.(subKey).crop_info = crop_info;

    catch ME
        warns{end+1} = sprintf('[WARN] %s: %s', subID, ME.message);
        Out.(subKey).data = [];
        Out.(subKey).crop_info = struct('error', ME.message);
    end
end

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
    fprintf('%s\n', strjoin(warns, newline));
end

end  % main function


%% ----------------- Local helpers -----------------

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
