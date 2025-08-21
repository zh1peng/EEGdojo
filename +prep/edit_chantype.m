function [EEG, out] = edit_chantype(EEG, varargin)
% EDIT_CHANTYPE Sets the type for each channel ('EEG', 'EOG', 'ECG', 'OTHER').
%
% This function classifies channels based on provided labels. Any channel not
% specified as EOG, ECG, or OTHER is defaulted to EEG. This is useful for
% constraining subsequent analyses (e.g., bad channel detection) to specific
% channel types, and for identifying artifactual components during ICA.
%
% Inputs:
%   EEG         - EEGLAB EEG structure.
%   varargin    - Optional parameters:
%     'EOGLabels'   - (cell array of strings) Labels of channels to be
%                     classified as Electrooculogram (EOG). Default is {}.
%     'ECGLabels'   - (cell array of strings) Labels of channels to be
%                     classified as Electrocardiogram (ECG). Default is {}.
%     'OtherLabels' - (cell array of strings) Labels of channels to be
%                     classified as 'OTHER'. Default is {}.
%     'LogFile'     - (char) Path to the log file for recording processing
%                     information. Default is ''.
%
% Outputs:
%   EEG         - Modified EEGLAB EEG structure with updated channel types.
%   out         - Structure containing output information:
%     .types_set  - (struct) A summary of how many channels were set to each type.
%
% Examples:
%   % 1. Classify EOG, ECG, and other channels:
%   EEG = edit_chantype(EEG, ...
%       'EOGLabels', {'VEOG', 'HEOG'}, ...
%       'ECGLabels', {'ECG1'}, ...
%       'OtherLabels', {'TRIG'});
%
%   % 2. Usage within a pipeline:
%   %    (Assuming 'p' is a parameter structure containing 'p.logFile')
%   pipe = pipe.addStep(@prep.edit_chantype, ...
%       'EOGLabels', {'VEOG', 'HEOG'}, ...
%       'ECGLabels', {'ECG1'}, ...
%       'OtherLabels', {'TRIG'}, ...
%       'LogFile', p.logFile);
%
% See also: chans2idx, eeg_checkset

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('EOGLabels', {}, @iscellstr);
    p.addParameter('ECGLabels', {}, @iscellstr);
    p.addParameter('OtherLabels', {}, @iscellstr);
    p.addParameter('LogFile', '', @(s) ischar(s) || isstring(s));

    p.parse(EEG, varargin{:});
    R = p.Results;

    out = struct('types_set', struct('EEG', 0, 'EOG', 0, 'ECG', 0, 'OTHER', 0));

    try
        logPrint(R.LogFile, '[edit_chantype] Starting channel type editing.');

        % Default all channels to 'EEG' first
        for i = 1:EEG.nbchan
            EEG.chanlocs(i).type = 'EEG';
        end
        
        % Set EOG channel types
        eog_idx = chans2idx(EEG, R.EOGLabels);
        for i = 1:length(eog_idx)
            EEG.chanlocs(eog_idx(i)).type = 'EOG';
        end
        if ~isempty(eog_idx)
            logPrint(R.LogFile, sprintf('[edit_chantype] Set %d channels to EOG: %s', length(eog_idx), strjoin(R.EOGLabels, ', ')));
        end

        % Set ECG channel types
        ecg_idx = chans2idx(EEG, R.ECGLabels);
        for i = 1:length(ecg_idx)
            EEG.chanlocs(ecg_idx(i)).type = 'ECG';
        end
        if ~isempty(ecg_idx)
            logPrint(R.LogFile, sprintf('[edit_chantype] Set %d channels to ECG: %s', length(ecg_idx), strjoin(R.ECGLabels, ', ')));
        end

        % Set OTHER channel types
        other_idx = chans2idx(EEG, R.OtherLabels);
        for i = 1:length(other_idx)
            EEG.chanlocs(other_idx(i)).type = 'OTHER';
        end
        if ~isempty(other_idx)
            logPrint(R.LogFile, sprintf('[edit_chantype] Set %d channels to OTHER: %s', length(other_idx), strjoin(R.OtherLabels, ', ')));
        end

        % Recalculate EEG channels (those not set to something else)
        all_non_eeg_idx = [eog_idx(:); ecg_idx(:); other_idx(:)]';
        eeg_idx = setdiff(1:EEG.nbchan, all_non_eeg_idx);

        % Summarize and log the final counts
        out.types_set.EOG = length(eog_idx);
        out.types_set.ECG = length(ecg_idx);
        out.types_set.OTHER = length(other_idx);
        out.types_set.EEG = length(eeg_idx);

        logPrint(R.LogFile, sprintf('[edit_chantype] Total channels classified: EEG=%d, EOG=%d, ECG=%d, OTHER=%d', ...
            out.types_set.EEG, out.types_set.EOG, out.types_set.ECG, out.types_set.OTHER));
        logPrint(R.LogFile, '[edit_chantype] Channel type editing complete.');
        
        EEG = eeg_checkset(EEG);

    catch ME
        error('[edit_chantype] Channel type editing failed: %s', ME.message);
    end
end