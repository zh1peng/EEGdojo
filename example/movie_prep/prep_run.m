%% Setup
EEGLAB_path = '/media/NAS/misc/matlab_toolbox/eeglab2023.1';
EEGdojo= '/media/NAS/misc/matlab_toolbox/EEGdojo';
FASTER_path =  '/media/NAS/misc/matlab_toolbox/FASTER';
addpath(genpath(EEGLAB_path))
addpath(genpath(FASTER_path))
addpath(genpath(EEGdojo))

clear; close all; clc;
%% Get subinfo
BIDS_path= '/media/NAS/EEGdata/BEAT';
task2prep = 'moviep';
derivative_path = sprintf('%s/derivatives/prep_%s',BIDS_path,task2prep);
if ~isfolder(derivative_path),mkdir(derivative_path),end;
subinfo_file =  fullfile(BIDS_path,'conversion_summary.tsv');
T = readtable(subinfo_file , 'FileType','text', 'Delimiter','\t');
param_file='/media/NAS/EEGdata/BEAT/code/Params.mat';

% Use output_set to derive everything
p = string(T.output_set);

% subid: supports sub-101 or sub_101
subid = regexp(p, '(sub[-_]\d+)', 'match', 'once');
% split into path and filename
[folders, names, exts] = arrayfun(@(s) fileparts(char(s)), p, ...
                                  'UniformOutput', false);
all_subpath  = folders;                   % e.g., .../sub-101/eeg
all_subfile  = strcat(names, exts);       % e.g., sub-101_task-rest_run-01_eeg.set

% extract task (e.g. rest, mid, sst, emodot, moviep…)
task_tokens = regexp(names, 'task-([^-_]+)', 'tokens', 'once');
task = cellfun(@(c) c{1}, task_tokens, 'UniformOutput', false);

% Example: filter for moviep
mask = strcmpi(task, task2prep);
subpath = all_subpath(mask);
subfile = all_subfile(mask);


%% Run preprocessing
tmp=load(param_file);

% parpool();
% check parallel pool
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool();

end
tic
parfor sub_i=1:length(subfile)
% Create a local copy of params
p=tmp.Params;
p.Input.filepath = subpath{sub_i};
p.Input.filename = subfile{sub_i};
p.Output.filepath = derivative_path;
p=setup_params(p);
logPrint(p.LogFile,'Starting ECG processing pipeline...');

% Define the ECG pipeline using the struct-passing convention
p.Output.filename = [regexprep(p.basename, 'eeg', 'ecg'), '.set'];
ecg_pipe = Pipeline([], p) ...
    .addStep('Load', @prep.load_set, p.Input) ...
    .addStep('SelectECG', @prep.select_channels, 'ChanLabels', p.ChanInfo.ECGChanLabel) ...
    .addStep('Downsample', @prep.downsample, p.Downsample) ...
    .addStep('Bandpass', @prep.filter, p.Filter) ...
    .addStep('Notch',@prep.remove_powerline, p.Powerline) ...
    .addStep('Crop by markers', @prep.crop_by_markers, p.Crop) ...
    .addStep('Save ECG', @prep.save_set, p.Output) ... % Override outfilename
    .run()

% Add the processed ECG data to the main parameter struct for the next pipeline
p.BadIC.ECG_Struct = ecg_pipe.EEG;
p.BadIC.DetectECG = true;
p.Output.filename = [p.basename, '.set'];
p.ChanInfo.Chan2remove = 'ECG';
% Define the main EEG pipeline using the struct-passing convention
eeg_pipe = prep.Pipeline([], p) ...
    .addStep('Load', @prep.load_set, p.Input) ...
    .addStep('Remove Other chans', @prep.remove_channels, 'Chan2remove', p.ChanInfo.Chan2remove) ...
    .addStep('Downsample', @prep.downsample, p.Downsample) ...
    .addStep('Bandpass', @prep.filter, p.Filter) ...
    .addStep('Notch',@prep.remove_powerline, p.Powerline) ...
    .addStep('Crop by markers', @prep.crop_by_markers, p.Crop) ...
    .addStep('Bad channels', @prep.remove_bad_channels, p.BadChan) ...
    .addStep('Re-reference', @prep.reref, p.Reref) ...
    .addStep('ICA + prune', @prep.remove_bad_ICs, p.BadIC) ...
    .addStep('Intepolate', @prep.interpolate, 'LogFile',p.LogFile) ...
    .addStep('Save', @prep.save_set, p.Output) ...
    .run()
logPrint(p.LogFile,'--- All pipelines finished ---');
end
toc