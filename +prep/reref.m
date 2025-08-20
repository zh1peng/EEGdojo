function EEG = reref(EEG, varargin)
%REREF Re-reference EEG data.
%   Usage: EEG = reref(EEG, 'excludeLabels', {'EOG1', 'EOG2'});

    p = inputParser;
    addRequired(p, 'EEG', @isstruct);
    addParameter(p, 'excludeLabels', {}, @(x) iscell(x) || ischar(x));
    addParameter(p, 'logFile', '', @ischar);
    parse(p, EEG, varargin{:});

    R = p.Results;
    EEG = R.EEG;

    if ischar(R.excludeLabels)
        excludeLabels = {R.excludeLabels};
    else
        excludeLabels = R.excludeLabels;
    end


    log_msg = '--- Re-referencing data to average';
    if ~isempty(excludeLabels)
        log_msg = [log_msg, sprintf(', excluding channels: %s', strjoin(excludeLabels, ', '))];
    end
    logPrint(R.logFile, [log_msg, ' ---']);


    if isempty(excludeLabels)
        EEG = pop_reref(EEG, []);
    else
        labels = {EEG.chanlocs.labels};
        excludeIdx = find(ismember(labels, excludeLabels));
        EEG = pop_reref(EEG, [], 'exclude', excludeIdx);
    end

    EEG = eeg_checkset(EEG);

    logPrint(R.logFile, 'Re-referencing complete.');

end
