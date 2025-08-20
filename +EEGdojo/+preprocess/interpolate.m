function EEG = interpolate(EEG, varargin)
%INTERPOLATE Interpolate missing channels.
%   Usage: EEG = interpolate(EEG);

    p = inputParser;
    addRequired(p, 'EEG', @isstruct);
    addParameter(p, 'logFile', '', @ischar);
    parse(p, EEG, varargin{:});

    R = p.Results;
    EEG = R.EEG;

    if ~isfield(EEG, 'urchanlocs') || isempty(EEG.urchanlocs)
        if ~isempty(R.logFile)
            logPrint(R.logFile, '--- Skipping interpolation: No original channel locations found. ---');
        end
        return;
    end

    original_chans = {EEG.urchanlocs.labels};
    current_chans = {EEG.chanlocs.labels};
    chans_to_interp = setdiff(original_chans, current_chans);

    if isempty(chans_to_interp)
        if ~isempty(R.logFile)
            logPrint(R.logFile, '--- Skipping interpolation: No channels to interpolate. ---');
        end
        return;
    end

    if ~isempty(R.logFile)
        logPrint(R.logFile, sprintf('--- Interpolating %d channels: %s ---', numel(chans_to_interp), strjoin(chans_to_interp, ', ')));
    end

    EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');
    EEG = eeg_checkset(EEG);

    if ~isempty(R.logFile)
        logPrint(R.logFile, 'Interpolation complete.');
    end
end
