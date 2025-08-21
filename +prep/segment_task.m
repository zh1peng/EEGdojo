function [EEG, out] = segment_task(EEG, varargin)
% SEGMENT_TASK  Segments continuous EEG data into epochs around specified event markers.
%   This function is designed for task-related EEG data where epochs are
%   created relative to specific time-locked event markers (e.g., stimulus
%   onset, response, trial start). It extracts segments of data within a
%   defined time window around these markers, preparing the data for
%   event-related potential (ERP) or other epoch-based analyses.
%
% Syntax:
%   [EEG, out] = prep.segment_task(EEG, 'param', value, ...)
%
% Input Arguments:
%   EEG         - EEGLAB EEG structure (continuous data with events).
%
% Optional Parameters (Name-Value Pairs):
%   'Markers'       - (cell array of strings, default: {})
%                     A cell array of event marker strings (e.g., {'S1', 'S2', 'Response'})
%                     around which to create epochs.
%   'TimeWindow'    - (numeric array [start_time end_time], default: [])
%                     A two-element numeric array specifying the time window
%                     for epoch extraction, relative to the event marker, in seconds.
%                     E.g., [-0.5 1.5] means from 500 ms before to 1500 ms after the marker.
%   'LogFile'       - (char | string, default: '')
%                     Path to a log file for verbose output. If empty, output
%                     is directed to the command window.
%
% Output Arguments:
%   EEG         - Modified EEGLAB EEG structure with data segmented into epochs.
%   out         - Structure containing details of the segmentation:
%                 out.epochs_created: A structure where field names are marker
%                                     types and values are the number of epochs
%                                     created for that marker.
%                 out.total_epochs: The total number of epochs created across
%                                   all specified markers.
%
% Examples:
%   % Example 1: Segment EEG around 'stim_on' and 'resp' markers (without pipeline)
%   % Load a continuous EEG dataset with events, e.g., EEG = pop_loadset('task_eeg.set');
%   [EEG_epoched, seg_info] = prep.segment_task(EEG, ...
%       'Markers', {'stim_on', 'resp'}, ...
%       'TimeWindow', [-0.2 0.8], ...
%       'LogFile', 'task_segmentation_log.txt');
%   disp('Task data segmented.');
%   disp('Epochs created per marker:');
%   disp(seg_info.epochs_created);
%   disp(['Total epochs: ', num2str(seg_info.total_epochs)]);
%
%   % Example 2: Segment with a different time window (with pipeline)
%   % Assuming 'pipe' is an initialized pipeline object
%   pipe = pipe.addStep(@prep.segment_task, ...
%       'Markers', {'trial_start'}, ...
%       'TimeWindow', [-1 2], ...
%       'LogFile', p.logFile); %% p.logFile from pipeline parameters
%   % Then run the pipeline: [EEG_processed, results] = pipe.run(EEG);
%   disp('Task data segmented via pipeline.');
%
% See also: pop_epoch

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('Markers', {}, @iscellstr);
    p.addParameter('TimeWindow', [], @(x) isnumeric(x) && numel(x) == 2);
    p.addParameter('LogFile', '', @(s) ischar(s) || isstring(s));

    p.parse(EEG, varargin{:});
    R = p.Results;

    out = struct();
    out.epochs_created = struct();

    if isempty(R.Markers) || isempty(R.TimeWindow)
        logPrint(R.LogFile, '[segment_task] Markers or TimeWindow is empty, skipping task segmentation.');
        return;
    end


    logPrint(R.LogFile, '[segment_task] ------ Segmenting task data ------');
    logPrint(R.LogFile, sprintf('[segment_task] Markers: %s, Time window: [%.2f %.2f]s', strjoin(R.Markers, ', '), R.TimeWindow(1), R.TimeWindow(2)));

    % Segment the data into epochs
    logPrint(R.LogFile, '[segment_task] Calling pop_epoch to segment data...');
    EEG = pop_epoch(EEG, R.Markers, R.TimeWindow, 'epochinfo', 'yes');
    
    % Log the number of epochs for each marker type
    unique_markers = unique(R.Markers);
    for i = 1:length(unique_markers)
        marker = unique_markers{i};
        n_epochs = sum(ismember({EEG.epoch.eventtype}, marker));
        out.epochs_created.(marker) = n_epochs;
        logPrint(R.LogFile, sprintf('[segment_task] Created %d epochs for marker %s', n_epochs, marker));
    end
    logPrint(R.LogFile, sprintf('[segment_task] Total epochs created: %d', EEG.trials));
    
    out.total_epochs = EEG.trials;

    logPrint(R.LogFile, '[segment_task] ------ Task segmentation complete ------');

end
