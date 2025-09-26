%% Setup
EEGLAB_path = '/media/NAS/misc/matlab_toolbox/eeglab2023.1';
EEGdojo     = '/media/NAS/misc/matlab_toolbox/EEGdojo';
FASTER_path = '/media/NAS/misc/matlab_toolbox/FASTER';
addpath(genpath(EEGLAB_path))
addpath(genpath(FASTER_path))
addpath(genpath(EEGdojo))

clear; close all; clc;

%% Get subinfo
BIDS_path   = '/media/NAS/EEGdata/BEAT';
tasks       = {'emodot'};   % <-- loop these
subinfo_file = fullfile(BIDS_path,'conversion_summary.tsv');
T = readtable(subinfo_file , 'FileType','text', 'Delimiter','\t');
param_file  = '/media/NAS/EEGdata/BEAT/code/emodot.mat';

% Use output_set to derive everything
p_all = string(T.output_set);

% subid: supports sub-101 or sub_101
subid = regexp(p_all, '(sub[-_]\d+)', 'match', 'once');

% split into path and filename
[all_paths, all_names, all_exts] = arrayfun(@(s) fileparts(char(s)), p_all, 'UniformOutput', false);
all_subpath  = all_paths;                         % e.g., .../sub-101/eeg
all_subfile  = strcat(all_names, all_exts);       % e.g., sub-101_task-rest_run-01_eeg.set

% extract task (e.g. rest, mid, sst, emodot, moviep…)
task_tokens = regexp(all_names, 'task-([^-_]+)', 'tokens', 'once');
all_task = cellfun(@(c) c{1}, task_tokens, 'UniformOutput', false);

%% Parallel pool (reuse across tasks)
poolobj = gcp('nocreate');
if isempty(poolobj), parpool(); end

%% Load template params once
tmp = load(param_file);   % expects tmp.Params

overallTic = tic;
for t = 1:numel(tasks)
    task2prep = tasks{t};
    fprintf('\n=== Processing task: %s ===\n', task2prep);

    % Per-task derivatives folder
    derivative_path = sprintf('%s/derivatives/prep_%s', BIDS_path, task2prep);
    if ~isfolder(derivative_path), mkdir(derivative_path); end

    % Filter subjects for this task
    mask    = strcmpi(all_task, task2prep);
    subpath = all_subpath(mask);
    subfile = all_subfile(mask);

    if isempty(subfile)
        fprintf('No entries found for task=%s. Skipping.\n', task2prep);
        continue;
    end

    

    % ---- Optional: per-task param tweaks (uncomment and adjust) ----
    % taskParams = tmp.Params;
    % switch lower(task2prep)
    %     case {'moviep','moviedm'}
    %         taskParams.Crop.StartMarker = '1';
    %         taskParams.Crop.EndMarker   = '2';
    %     case 'sst'
    %         % taskParams.Crop = ...
    %     case 'mid'
    %         % taskParams.Crop = ...
    % end

    tic
    parfor i = 1:numel(subfile)
        % Create a local copy of params
        p = tmp.Params;               % or taskParams if you enable the block above
        p.Input.filepath  = subpath{i};
        p.Input.filename  = subfile{i};
        p.Output.filepath = derivative_path;

        % Setup logging + derived names
        p = setup_params(p);
        logPrint(p.LogFile, sprintf('Starting pipelines for %s | %s', task2prep, subfile{i}));

        % --------- ECG pipeline ---------
        p_ecg = p;
        p_ecg.Output.filename = [regexprep(p_ecg.basename, 'eeg', 'ecg'), '.set'];

        ecg_pipe = prep.Pipeline([], p_ecg) ...
            .addStep('Load',              @prep.load_set,          p_ecg.Input) ...
            .addStep('SelectECG',         @prep.select_channels, 'ChanLabels', p_ecg.ChanInfo.ECGChanLabel) ...
            .addStep('Downsample',        @prep.downsample,      p_ecg.Downsample) ...
            .addStep('Bandpass',          @prep.filter,          p_ecg.Filter) ...
            .addStep('Notch',             @prep.remove_powerline,p_ecg.Powerline) ...
            .addStep('Crop by markers',   @prep.crop_by_markers, p_ecg.Crop) ...
            .addStep('Save ECG',          @prep.save_set,          p_ecg.Output) ...
            .run();

        % --------- EEG pipeline ---------
        p_eeg = p;
        p_eeg.BadIC.ECG_Struct   = ecg_pipe.EEG;
        p_eeg.BadIC.DetectECG    = true;
        p_eeg.Output.filename    = [p_eeg.basename, '.set'];
        p_eeg.ChanInfo.Chan2remove = 'ECG';  % drop ECG from EEG

        eeg_pipe = prep.Pipeline([], p_eeg) ...
            .addStep('Load',               @prep.load_set,          p_eeg.Input) ...
            .addStep('Remove Other chans', @prep.remove_channels, 'Chan2remove', p_eeg.ChanInfo.Chan2remove) ...
            .addStep('Downsample',         @prep.downsample,      p_eeg.Downsample) ...
            .addStep('Bandpass',           @prep.filter,          p_eeg.Filter) ...
            .addStep('Notch',              @prep.remove_powerline,p_eeg.Powerline) ...
            .addStep('Crop by markers',    @prep.crop_by_markers, p_eeg.Crop) ...
            .addStep('Bad channels',       @prep.remove_bad_channels, p_eeg.BadChan) ...
            .addStep('Re-reference',       @prep.reref,           p_eeg.Reref) ...
            .addStep('ICA + prune',        @prep.remove_bad_ICs,  p_eeg.BadIC) ...
            .addStep('Intepolate',         @prep.interpolate,     'LogFile', p_eeg.LogFile) ...
            .addStep('Save',               @prep.save_set,          p_eeg.Output) ...
            .run();

        logPrint(p.LogFile, sprintf('Finished %s | %s', task2prep, subfile{i}));
    end
    fprintf('Task %s done. Elapsed: %.1f s\n', task2prep, toc);
end
fprintf('\nAll tasks finished. Total elapsed: %.1f s\n', toc(overallTic));
