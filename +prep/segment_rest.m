function [EEG, out] = segment_rest(EEG, varargin)
% SEGMENT_REST  Segments continuous EEG data into overlapping epochs for resting-state analysis.
%   This function is specifically designed for continuous EEG data that does
%   not have pre-defined time-locked events, such as resting-state recordings.
%   It works by first generating regularly spaced 'epoch_start' events within
%   the continuous data. These synthetic events are then used to create
%   overlapping epochs of a specified length. This method is useful for
%   preparing continuous resting-state data for subsequent epoch-based analyses.
%
% Syntax:
%   [EEG, out] = prep.segment_rest(EEG, 'param', value, ...)
%
% Input Arguments:
%   EEG         - EEGLAB EEG structure (continuous data).
%
% Optional Parameters (Name-Value Pairs):
%   'EpochLength'   - (numeric, default: 2)
%                     The desired length of each epoch in seconds.
%   'EpochOverlap'  - (numeric, default: 0.5)
%                     The desired overlap between consecutive epochs, as a
%                     proportion (0 to <1). For example, 0.5 means 50% overlap.
%                     A value of 0 means no overlap (contiguous epochs).
%   'LogFile'       - (char | string, default: '')
%                     Path to a log file for verbose output. If empty, output
%                     is directed to the command window.
%
% Output Arguments:
%   EEG         - Modified EEGLAB EEG structure with data segmented into epochs.
%   out         - Structure containing details of the segmentation:
%                 out.epochs_created: The number of epochs created.
%
% Examples:
%   % Example 1: Segment continuous EEG into 2-second epochs with 50% overlap (without pipeline)
%   % Load a continuous EEG dataset first, e.g., EEG = pop_loadset('continuous_eeg.set');
%   [EEG_epoched, seg_info] = prep.segment_rest(EEG, ...
%       'EpochLength', 2, ...
%       'EpochOverlap', 0.5, ...
%       'LogFile', 'segmentation_log.txt');
%   disp(['Created ', num2str(seg_info.epochs_created), ' epochs.']);
%
%   % Example 2: Segment into 4-second epochs with 25% overlap (with pipeline)
%   % Assuming 'pipe' is an initialized pipeline object
%   pipe = pipe.addStep(@prep.segment_rest, ...
%       'EpochLength', 4, ...
%       'EpochOverlap', 0.25, ...
%       'LogFile', p.logFile); %% p.logFile from pipeline parameters
%   % Then run the pipeline: [EEG_processed, results] = pipe.run(EEG);
%   disp('Resting-state data segmented via pipeline.');
%
% See also: eeg_regepochs, pop_epoch

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('EpochLength', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('EpochOverlap', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
    p.addParameter('LogFile', '', @(s) ischar(s) || isstring(s));

    p.parse(EEG, varargin{:});
    R = p.Results;

    out = struct();

    if isempty(R.EpochLength)
        logPrint(R.LogFile, '[segment_rest] EpochLength is empty, skipping segmentation.');
        return;
    end


    logPrint(R.LogFile, '[segment_rest] ------ Segmenting resting-state data ------');
    logPrint(R.LogFile, sprintf('[segment_rest] Epoch length: %.2f s, Overlap: %.2f%%', R.EpochLength, R.EpochOverlap*100));

    % Create regularly spaced markers for epoching
    logPrint(R.LogFile, '[segment_rest] Creating regularly spaced "epoch_start" events...');
    EEG = eeg_regepochs(EEG, 'recurrence', R.EpochLength * (1 - R.EpochOverlap), 'eventtype', 'epoch_start', 'extractepochs', 'off');
    logPrint(R.LogFile, sprintf('[segment_rest] %d "epoch_start" events created.', length(EEG.event)));

    % Segment the data into epochs based on the new markers
    logPrint(R.LogFile, '[segment_rest] Epoching data based on "epoch_start" events...');
    EEG = pop_epoch(EEG, {'epoch_start'}, [0 R.EpochLength], 'epochinfo', 'yes');

    out.epochs_created = EEG.trials;

    logPrint(R.LogFile, sprintf('[segment_rest] Data segmented into %d epochs.', out.epochs_created));
    logPrint(R.LogFile, '[segment_rest] ------ Segmentation complete ------');

end
