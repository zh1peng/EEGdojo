function Out = extract_segment(varargin)
% EXTRACT_SEGMENT Build a study-level container from cropped segments of many EEGLAB .set files.
%   This function is designed to find EEGLAB .set files and extract a single,
%   continuous data segment from each file based on specified start and end
%   markers. It's ideal for analyzing data from continuous recordings like
%   movie-watching or resting-state paradigms.
%
% Output Structure:
%   The output is a single struct `Out` with two main components:
%   1. Out.meta: A struct containing shared information across all subjects.
%      This includes:
%       .fs              - Sampling rate of the processed data.
%       .srate           - (Same as fs) Sampling rate.
%       .chanlocs        - Channel locations structure from the first subject.
%       .start_marker    - The event marker used to define the start of the segment.
%       .end_marker      - The event marker used to define the end of the segment.
%       .padding_sec     - Seconds of padding added before the start and after the end marker.
%       .created_at      - Timestamp of when the extraction was run.
%       .run_tag         - A user-defined tag for this specific extraction run.
%
%   2. Out.sub_...: A series of dynamically named fields for each subject
%      (e.g., Out.sub_001, Out.sub_S02). Each subject field is a struct containing:
%       .data            - A [channels x time] matrix of the extracted segment.
%       .crop_info       - A struct with info from `prep.crop_by_markers`,
%                          including the start and end sample points of the crop.
%
% Name-Value Pair Arguments:
%   'study_path' (char, REQUIRED)
%       The absolute path to the root directory containing the .set files.
%
%   'StartMarker' (char/string, REQUIRED)
%       The event marker string that defines the beginning of the segment to extract.
%
%   'EndMarker' (char/string, REQUIRED)
%       The event marker string that defines the end of the segment.
%
%   'searchstring' (char, default: '.*\.set$')
%       A regular expression to find specific .set files within the study_path.
%
%   'recursive' (logical, default: true)
%       If true, searches for files in subdirectories of study_path.
%
%   'subject_parser' (char, default: '(?<sub>.+)')
%       A regular expression with a named token '(?<sub>...)' to extract a unique
%       subject ID from each filename.
%
%   'chan_include' (cellstr/string, default: {})
%       A list of channel labels to *include*. All other channels will be removed.
%
%   'chan_exclude' (cellstr/string, default: {})
%       A list of channel labels to *exclude*. All other channels will be kept.
%
%   'PadSec' (numeric, default: 0)
%       Number of seconds to add as padding before the StartMarker and after the EndMarker.
%
%   'to_single' (logical, default: false)
%       If true, casts the output data matrices to the 'single' data type to save memory.
%
%   'save_path' (char, default: '')
%       If provided, the final 'Out' struct is saved to this full path.
%
%   'save_v7_3' (logical, default: true)
%       If saving, use the '-v7.3' flag to support larger files.
%
%   'run_tag' (char, default: timestamp)
%       A custom string tag to identify this extraction run, stored in Out.meta.
%
% Examples:
%
%   % Example 1: Basic Extraction
%   % Extracts the segment between 'movie_start' and 'movie_end' from all
%   % .set files in the specified directory.
%   MovieData = Movie.extract_segment( ...
%       'study_path', 'C:\EEG_Data\MovieStudy',
%       'StartMarker', 'movie_start',
%       'EndMarker', 'movie_end');
%
%   % Example 2: Resting-State with Padding and Saving
%   % Extracts a resting-state segment, adds 5 seconds of padding to each
%   % side, and saves the output to a file.
%   RestData = Movie.extract_segment( ...
%       'study_path', 'C:\EEG_Data\RestingState',
%       'StartMarker', 'start_rest',
%       'EndMarker', 'end_rest',
%       'PadSec', 5,
%       'save_path', 'C:\EEG_Data\RestingState\processed\rest_segments.mat',
%       'run_tag', 'resting_state_v1');
%
%   % Example 3: Excluding Noisy Channels
%   % Extracts a segment but removes EOG and other specified channels first.
%   CleanData = Movie.extract_segment( ...
%       'study_path', 'C:\EEG_Data\MovieStudy',
%       'StartMarker', 'start',
%       'EndMarker', 'stop',
%       'chan_exclude', {'HEOG', 'VEOG', 'M1', 'M2'});

%% ---- Parse Name–Value inputs
p = inputParser;
p.FunctionName = 'extract_segment';
addParameter(p, 'study_path', '', @(x) ischar(x) && ~isempty(x));
addParameter(p, 'searchstring', '.*\.set$', @ischar);
addParameter(p, 'recursive', true, @islogical);
addParameter(p, 'subject_parser', '(?<sub>.+)', @ischar);

addParameter(p,'chan_include',{},@(x)iscellstr(x) || isstring(x));
addParameter(p,'chan_exclude',{},@(x)iscellstr(x) || isstring(x));

addParameter(p, 'StartMarker', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'EndMarker', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'PadSec', 0, @isnumeric);

addParameter(p, 'locutoff', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'hicutoff', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'csd', false, @islogical);
addParameter(p, 'csd_lambda', 1.0e-5, @isscalar);
addParameter(p, 'csd_head_radius', 10.0, @isscalar);
addParameter(p, 'csd_m', 4, @isscalar);
addParameter(p, 'hilbert_env', false, @islogical);

% Time-frequency decomposition
addParameter(p, 'tfr', false, @islogical);
addParameter(p, 'tfr_freqs', [], @(x) isempty(x) || isvector(x));
addParameter(p, 'tfr_n_cycles', [3 0.5], @isnumeric);

addParameter(p, 'to_single', false, @islogical);
addParameter(p, 'save_path', '', @ischar);
addParameter(p, 'save_v7_3', true, @islogical);
addParameter(p, 'run_tag', datestr(now, 'yyyy-mm-dd_HHMMSS'), @ischar);

parse(p, varargin{:});
opt = p.Results;

assert(~isempty(opt.study_path) && isfolder(opt.study_path), 'study_path not found.');
assert(~isempty(opt.StartMarker) && ~isempty(opt.EndMarker), 'StartMarker and EndMarker are required.');
if opt.tfr
    assert(exist('newtimef', 'file'), 'newtimef() not found. Is EEGLAB in your path?');
    assert(~isempty(opt.tfr_freqs), 'tfr_freqs must be provided when tfr=true.');
end


%% ---- Find .set files
[paths, names] = filesearch_regexp(opt.study_path, opt.searchstring, opt.recursive);
setFiles = fullfile(paths, names);
assert(~isempty(setFiles), 'No .set files matched searchstring in study_path.');

%% ---- Initialize output container
Out = struct();
preproc_opts = struct('locutoff', opt.locutoff, 'hicutoff', opt.hicutoff, ...
    'csd', opt.csd, 'csd_lambda', opt.csd_lambda, ...
    'csd_head_radius', opt.csd_head_radius, 'csd_m', opt.csd_m, ...
    'hilbert_env', opt.hilbert_env);

tfr_opts = struct('enable', opt.tfr, 'freqs', opt.tfr_freqs, 'n_cycles', opt.tfr_n_cycles);

Out.meta = struct( ...
    'fs',              [], ...
    'chanlocs',        [], ...
    'start_marker',    opt.StartMarker, ...
    'end_marker',      opt.EndMarker, ...
    'padding_sec',     opt.PadSec, ...
    'preproc',         preproc_opts, ...
    'tfr',             tfr_opts, ...
    'tfr_freqs',       [], ...
    'tfr_times',       [], ...
    'created_at',      datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
    'run_tag',         opt.run_tag, ...
    'version',         1.1);

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

        % ---- Channel selection
        if ~isempty(opt.chan_include) && ~isempty(opt.chan_exclude)
            error('Use either chan_include or chan_exclude, not both.');
        end
        if ~isempty(opt.chan_include)
            EEG = pop_select(EEG, 'channel', opt.chan_include);
            EEG = eeg_checkset(EEG);
        end
        if ~isempty(opt.chan_exclude)
            EEG = pop_select(EEG, 'nochannel', opt.chan_exclude);
            EEG = eeg_checkset(EEG);
        end

        % ---- Filtering
        if (~isempty(opt.locutoff) && opt.locutoff > 0) || (~isempty(opt.hicutoff) && opt.hicutoff > 0)
            args = {};
            if ~isempty(opt.locutoff) && opt.locutoff > 0, args = [args, {'locutoff', opt.locutoff}]; end
            if ~isempty(opt.hicutoff) && opt.hicutoff > 0, args = [args, {'hicutoff', opt.hicutoff}]; end
            EEG = pop_eegfiltnew(EEG, args{:});
            EEG = eeg_checkset(EEG);
        end

        % ---- CSD Transformation
        if opt.csd
            if ~exist('CSD.m', 'file')
                error('CSD.m not found. Make sure the CSD toolbox is in your MATLAB path.');
            end
            montage.lab = {EEG.chanlocs.labels}';
            montage.theta = [EEG.chanlocs.sph_theta]';
            montage.phi = [EEG.chanlocs.sph_phi]';
            [G, H] = GetGH(montage, opt.csd_m);
            EEG.data = CSD(EEG.data, G, H, opt.csd_lambda, opt.csd_head_radius);
        end

        % ---- Hilbert Amplitude Envelope
        if opt.hilbert_env
            for chan_i = 1:EEG.nbchan
                EEG.data(chan_i, :) = abs(hilbert(EEG.data(chan_i, :)));
            end
        end

        % Crop the segment
        [EEG_cropped, crop_info] = prep.crop_by_markers(EEG, ...
            'StartMarker', opt.StartMarker, ...
            'EndMarker', opt.EndMarker, ...
            'PadSec', opt.PadSec);

        % ---- Time-Frequency Decomposition
        if opt.tfr
            fprintf('Running TFR for %s...\n', subID);
            % Per user insight, newtimef must be called per-channel for continuous data
            % to avoid misinterpreting the [chan x pnts] matrix.
            
            % Run the first channel to get output dimensions
            [ersp_one_chan, ~, ~, tfr_times, tfr_freqs] = newtimef(EEG_cropped.data(1,:), ...
                EEG_cropped.pnts, ...
                [EEG_cropped.xmin EEG_cropped.xmax]*1000, ...
                EEG_cropped.srate, ...
                opt.tfr_n_cycles, ...
                'winsize', 250, ...
                'freqs', opt.tfr_freqs, ...
                'baseline', NaN, ...
                'plotersp', 'off', ...
                'plotitc', 'off', ...
                'verbose', 'off');

            % Pre-allocate the full data matrix
            data = zeros(EEG_cropped.nbchan, size(ersp_one_chan,1), size(ersp_one_chan,2));
            data(1,:,:) = ersp_one_chan;

            % Loop through the rest of the channels
            for chan_idx = 2:EEG_cropped.nbchan
                [data(chan_idx,:,:)] = newtimef(EEG_cropped.data(chan_idx,:), ...
                    EEG_cropped.pnts, ...
                    [EEG_cropped.xmin EEG_cropped.xmax]*1000, ...
                    EEG_cropped.srate, ...
                    opt.tfr_n_cycles, ...
                     'winsize', 250, ...
                    'freqs', opt.tfr_freqs, ...
                    'baseline', NaN, ...
                    'plotersp', 'off', ...
                    'plotitc', 'off', ...
                    'verbose', 'off');
            end
        else
            data = EEG_cropped.data;
        end

        % Lock meta from the first processed file
        if ~firstMetaLocked
            Out.meta.fs = EEG_cropped.srate;
            Out.meta.chanlocs = EEG_cropped.chanlocs;
            if opt.tfr
                Out.meta.tfr_freqs = tfr_freqs;
                Out.meta.tfr_times = tfr_times;
            end
            firstMetaLocked = true;
        else
            if EEG_cropped.srate ~= Out.meta.fs
                warns{end+1} = sprintf('[WARN] %s: srate %g ~= meta.srate %g', subID, EEG_cropped.srate, Out.meta.fs);
            end
            if numel(EEG_cropped.chanlocs) ~= numel(Out.meta.chanlocs)
                warns{end+1} = sprintf('[WARN] %s: channel count differs from first file', subID);
            end
        end

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