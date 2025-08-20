function EEG = filter(EEG, varargin)
%FILTER Apply high-pass and low-pass filters to EEG data.
%   Usage: EEG = filter(EEG, 'HPfreq', 0.5, 'LPfreq', 40);

    p = inputParser;
    addRequired(p, 'EEG', @isstruct);
    addParameter(p, 'HPfreq', -1, @isnumeric);
    addParameter(p, 'LPfreq', -1, @isnumeric);
    addParameter(p, 'logFile', '', @ischar);
    parse(p, EEG, varargin{:});

    R = p.Results;
    EEG = R.EEG;


    log_msg = '--- Applying filter --- \n ';
    if R.HPfreq > 0, log_msg = [log_msg, sprintf('High-pass: %.2f Hz. ', R.HPfreq)]; end
    if R.LPfreq > 0, log_msg = [log_msg, sprintf('Low-pass: %.2f Hz.', R.LPfreq)]; end
    logPrint(R.logFile, log_msg);
    
    if R.HPfreq < R.LPfreq
        error('[filter] High-pass frequency must be lower than low-pass.');
    end

    if R.HPfreq < 0 || R.LPfreq < 0
        error('[filter] Frequencies must be positive.');
    end

    try
        EEG = pop_eegfiltnew(EEG, 'locutoff', R.HPfreq, 'plotfreqz', 0);
        EEG = pop_eegfiltnew(EEG, 'hicutoff', R.LPfreq,  'plotfreqz', 0);
        EEG = eeg_checkset(EEG);
         logPrint(R.logFile, '--- Filter complete. ---');
    catch ME
        error('[filter] Filter failed: %s.',ME.message);
    end 
end
