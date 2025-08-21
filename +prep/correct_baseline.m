function [EEG, out] = correct_baseline(EEG, varargin)
% CORRECT_BASELINE Performs baseline correction on EEG data.
%
% This function applies baseline correction to EEG data, either epoched
% or continuous, by subtracting the mean of a specified baseline window.
% It uses the EEGLAB function `pop_rmbase`.
%
% Inputs:
%   EEG         - EEGLAB EEG structure.
%   varargin    - Optional parameters:
%     'BaselineWindow' - (numeric array) A two-element array [start_ms, end_ms]
%                        specifying the baseline window in milliseconds.
%                        e.g., [-200 0] for 200ms before stimulus onset.
%                        Default is []. If empty, baseline correction is skipped.
%     'LogFile'        - (char) Path to the log file for recording processing
%                        information. Default is ''.
%
% Outputs:
%   EEG         - Modified EEGLAB EEG structure with baseline corrected data.
%   out         - Structure containing output information:
%     .baseline_window_ms - (numeric array) The baseline window used for correction.
%
% Examples:
%   % 1. Apply baseline correction from -200ms to 0ms:
%   EEG = correct_baseline(EEG, 'BaselineWindow', [-200 0]);
%
%   % 2. Usage within a pipeline:
%   %    (Assuming 'p' is a parameter structure containing 'p.logFile')
%   pipe = pipe.addStep(@prep.correct_baseline, ...
%       'BaselineWindow', [-500 0], ...
%       'LogFile', p.logFile);
%
% See also: pop_rmbase

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('BaselineWindow', [], @(x) isnumeric(x) && numel(x) == 2);
    p.addParameter('LogFile', '', @(s) ischar(s) || isstring(s));

    p.parse(EEG, varargin{:});
    R = p.Results;

    out = struct();
    out.baseline_window_ms = R.BaselineWindow;

    if isempty(R.BaselineWindow)
        logPrint(R.LogFile, '[correct_baseline] BaselineWindow is empty, skipping baseline correction.');
        return;
    end

    try
        logPrint(R.LogFile, '[correct_baseline] Starting baseline correction.');
        logPrint(R.LogFile, sprintf('[correct_baseline] Baseline window: [%d %d] ms', R.BaselineWindow(1), R.BaselineWindow(2)));

        if EEG.trials > 1
            logPrint(R.LogFile, '[correct_baseline] Applying baseline correction to epoched data.');
        else
            logPrint(R.LogFile, '[correct_baseline] Applying baseline correction to continuous data.');
        end

        % Perform baseline correction
        EEG = pop_rmbase(EEG, R.BaselineWindow);
        
        logPrint(R.LogFile, '[correct_baseline] Baseline correction complete.');

    catch ME
        error('[correct_baseline] Baseline correction failed: %s', ME.message);
    end
end