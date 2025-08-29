function [EEG, out] = remove_bad_epoch(EEG, varargin)
% REMOVE_BAD_EPOCH  Detects and removes bad epochs from an EEG dataset.
%   This function identifies bad epochs using EEGLAB's `pop_autorej` and
%   the FASTER algorithm's `epoch_properties`. It then removes the
%   identified bad epochs from the dataset.
%
% Syntax:
%   [EEG, out] = prep.remove_bad_epoch(EEG, 'param', value, ...)
%
% Input Arguments:
%   EEG         - EEGLAB EEG structure (epoched data).
%
% Optional Parameters (Name-Value Pairs):
%   'Autorej'           - (logical, default: true)
%                         Enable epoch rejection using `pop_autorej`.
%   'Autorej_MaxRej'    - (numeric, default: 2)
%                         `maxrej` parameter for `pop_autorej`.
%   'FASTER'            - (logical, default: true)
%                         Enable epoch rejection using FASTER's `epoch_properties`.
%   'LogFile'           - (char | string, default: '')
%                         File path to log the results.
%
% Output Arguments:
%   EEG         - Modified EEGLAB EEG structure with bad epochs removed.
%   out         - Structure containing details of the detection:
%                 out.Bad.autorej: Indices of bad epochs from `pop_autorej`.
%                 out.Bad.FASTER: Indices of bad epochs from FASTER.
%                 out.Bad.all: Combined unique indices of all bad epochs.
%                 out.summary: Summary of bad epochs per detector and total.
%
% See also: pop_autorej, epoch_properties, pop_select

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('Autorej', true, @islogical);
    p.addParameter('Autorej_MaxRej', 2, @isnumeric);
    p.addParameter('FASTER', true, @islogical);
    p.addParameter('LogFile', '', @(s) ischar(s) || isstring(s));
    p.parse(EEG, varargin{:});
    R = p.Results;

    out = struct(); % Initialize out struct

    logPrint(R.LogFile, 'Identifying bad epochs...');

    % ----------------- Run Detectors -----------------
    Bad = struct();
    
    if R.Autorej
        logPrint(R.LogFile, '[remove_bad_epoch] Running pop_autorej detector...');
        [~, Bad.autorej] = pop_autorej(EEG, 'nogui', 'on', 'maxrej', R.Autorej_MaxRej);
    else
        Bad.autorej = [];
    end

    if R.FASTER
        logPrint(R.LogFile, '[remove_bad_epoch] Running FASTER detector...');
        epoch_list = epoch_properties(EEG, 1:EEG.nbchan);
        Bad.FASTER = find(min_z(epoch_list) == 1)';
    else
        Bad.FASTER = [];
    end

    % ----------------- Combine & Summarize -----------------
    Bad.all = unique([Bad.autorej, Bad.FASTER]);

    summary = struct();
    summary.autorej = numel(Bad.autorej);
    summary.FASTER = numel(Bad.FASTER);
    summary.Total = numel(Bad.all);

    logPrint(R.LogFile, sprintf('Bad Epochs Identified by Auto Reject: %d\nDetails: %s', summary.autorej, mat2str(Bad.autorej)));
    logPrint(R.LogFile, sprintf('Bad Epochs Identified by FASTER: %d\nDetails: %s', summary.FASTER, mat2str(Bad.FASTER)));
    logPrint(R.LogFile, sprintf('Total Unique Bad Epochs: %d\nDetails: %s\n', summary.Total, mat2str(Bad.all)));

    % ----------------- Action: Remove -----------------
    if ~isempty(Bad.all)
        logPrint(R.LogFile, sprintf('[remove_bad_epoch] Removing %d bad epochs...', summary.Total));
        EEG = pop_select(EEG, 'notrial', Bad.all);
        EEG = eeg_checkset(EEG);
        logPrint(R.LogFile, '[remove_bad_epoch] Bad epochs removed successfully.');
    else
        logPrint(R.LogFile, '[remove_bad_epoch] No bad epochs to remove.');
    end

    % ----------------- Bookkeeping in EEG.etc -----------------
    out.Bad = Bad;
    out.summary = summary;

    if ~isfield(EEG.etc, 'EEGdojo'), EEG.etc.EEGdojo = struct(); end
    EEG.etc.EEGdojo.BadEpochIdx = Bad.all;
    EEG.etc.EEGdojo.BadEpochSummary = summary;

end