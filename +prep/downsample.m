function EEG = downsample(EEG, varargin)
%DOWNSAMPLE Downsample EEG data.
%   Usage: EEG = downsample(EEG, 'freq', 250);

    p = inputParser;
    addRequired(p, 'EEG', @isstruct);
    addParameter(p, 'freq', 250, @isnumeric);
    addParameter(p, 'logFile', '', @ischar);
    addParameter(p, 'error_logFile', '', @ischar);
    parse(p, EEG, varargin{:});

    
    R = p.Results;
    EEG = R.EEG;
    logPrint(R.logFile, sprintf('--- Downsampling data to %d Hz ---', R.freq));
    try
        EEG = pop_resample(EEG, R.freq);
        EEG = eeg_checkset(EEG);
        logPrint(R.logFile, sprintf('Downsampling complete. New sampling rate: %d Hz.', EEG.srate));
    catch ME 
        logPrint(R.error_logFile, 'Downsampling failed: %s.',ME.message);
        return;
    end
end
