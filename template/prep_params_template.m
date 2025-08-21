params=struct;

% General I/O and Logging parameters
params.IO.filename=''; % Example: 'sub-NDARAD481FXF_task-RestingState_eeg.set';
params.IO.filepath=''; % Example: 'V:\HBN_EEG_BIDS\cmi_bids_NC\sub-NDARAD481FXF\eeg';
params.IO.outputpath=''; % Example: 'X:\HBN_EEG_preprocessing\code\test';
% params.IO = setIOPaths(params.IO); % Uncomment and configure setIOPaths if needed
params.Log.LogPath = pwd; % Path to save log files and plots
params.Log.LogFile = '';  % Base name for the log report file (e.g., 'preprocessing_log')

% Channel information
params.ChanInfo.AllChan = [1:64];
params.ChanInfo.RefChan = 129; % Reference channel index (if applicable)
params.ChanInfo.EOGChanIdx = []; % Indices of EOG channels
params.ChanInfo.EOGChanLabels = {}; % Labels of EOG channels (e.g., {'VEOG', 'HEOG'})
params.ChanInfo.KnownBadChanIdx = []; % Indices of channels known to be bad
params.ChanInfo.KnownBadChanLabels = {}; % Labels of channels known to be bad
params.ChanInfo.Chan2remove = {}; % Channels to remove (e.g., {'Cz', 'Fz'} or [1 5 10])

% Downsampling and filtering parameters
params.DownsamplingRate = 250; % Target downsampling rate in Hz
params.Filter.LowCutoff = 0.1; % High-pass filter cutoff in Hz
params.Filter.HighCutoff = 40; % Low-pass filter cutoff in Hz

% Powerline noise removal parameters (prep.remove_powerline)
params.Powerline.Method = 'cleanline'; % 'cleanline' or 'notch'
params.Powerline.Freq = 60; % Fundamental powerline frequency (e.g., 50 or 60 Hz)
params.Powerline.BW = 2; % Half-bandwidth for FIR notch filter (only for 'notch' method)
params.Powerline.NHarm = 3; % Number of harmonics to target

% Bad channel detection parameters (prep.remove_bad_channels)
params.BadChan.IdxDetect = setdiff(params.ChanInfo.AllChan, [params.ChanInfo.RefChan, params.ChanInfo.EOGChanIdx, params.ChanInfo.KnownBadChanIdx]); % Channels to consider for detection
params.BadChan.Action = 'remove'; % 'remove' or 'flag'
params.BadChan.KnownBadIdx = params.ChanInfo.KnownBadChanIdx; % Indices of channels already known to be bad

% Classic EEGLAB detectors (pop_rejchan)
params.BadChan.Kurtosis = false; % Enable kurtosis-based detection
params.BadChan.Kurt_Threshold = 5; % Kurtosis threshold (z-score)
params.BadChan.Probability = false; % Enable probability-based detection
params.BadChan.Prob_Threshold = 5; % Probability threshold (z-score)
params.BadChan.Spectrum = false; % Enable spectrum-based detection
params.BadChan.Spec_Threshold = 5; % Spectrum threshold (z-score)
params.BadChan.Spec_FreqRange = [0.1 40]; % Frequency range for spectrum analysis (Hz)
params.BadChan.NormOn = 'on'; % Normalize measures ('on' or 'off')

% FASTER detectors (FASTER_rejchan)
params.BadChan.FASTER_MeanCorr = false; % Enable FASTER mean correlation
params.BadChan.FASTER_Threshold = 0.4; % FASTER mean correlation threshold
params.BadChan.FASTER_Variance = false; % Enable FASTER variance
params.BadChan.FASTER_VarThreshold = 3; % FASTER variance threshold
params.BadChan.FASTER_Hurst = false; % Enable FASTER Hurst exponent
params.BadChan.FASTER_HurstThreshold = 3; % FASTER Hurst exponent threshold

% CleanRaw detectors (cleanraw_rejchan)
params.BadChan.CleanRaw_Flatline = true; % Enable CleanRaw flatline detection
params.BadChan.Flatline_Sec = 5; % Minimum flatline duration in seconds
params.BadChan.CleanDrift_Band = [0.25 0.75]; % Highpass filter band for CleanRaw methods
params.BadChan.CleanRaw_Noise = true; % Enable CleanRaw noise detection (channel correlation)
params.BadChan.CleanChan_Corr = 0.8; % Correlation threshold for CleanRaw noise
params.BadChan.CleanChan_Line = 4; % Line noise threshold for CleanRaw noise
params.BadChan.CleanChan_MaxBad = 0.5; % Maximum proportion of bad time points for CleanRaw noise
params.BadChan.CleanChan_NSamp = 50; % Number of samples for CleanRaw noise

% Re-referencing parameters (prep.reref)
params.Reref.excludeLabels = setdiff(params.ChanInfo.AllChan, [params.ChanInfo.RefChan, params.ChanInfo.EOGChanIdx, params.ChanInfo.KnownBadChanIdx]); % Channels to exclude from average reference
% params.Reref.Chan2averageOn = 'Off'; % Custom flag for pipeline if needed

% ICA parameters (prep.remove_bad_ICs)
params.BadIC.FilterICAOn = true; % Apply high-pass filter before ICA
params.BadIC.FilterICALocutoff = 1; % High-pass filter cutoff for ICA (Hz)
params.BadIC.ICAType = 'runica'; % Type of ICA algorithm
params.BadIC.ICLabelOn = true; % Enable ICLabel classification
params.BadIC.ICLabelThreshold = [0 0.1; 0.9 1]; % ICLabel thresholds
params.BadIC.FASTEROn = true; % Enable FASTER component property analysis
params.BadIC.DetectECG = false; % Enable ECG correlation-based detection
params.BadIC.ECGStruct = []; % New: Placeholder for external ECG data structure
params.BadIC.ECGCorrelationThreshold = 0.8; % ECG correlation threshold

% Bad point detection parameters (clean_rawdata) - Not directly modified by prep functions
params.BadPoint.on = "on";
params.BadPoint.BurstCriterion = 20;
params.BadPoint.WindowCriterion = 0.25;
params.BadPoint.WindowCriterionTolerances = [-Inf 7];

% Epoching resting state data parameters (prep.segment_rest)
params.EpochRest.EpochLength = 2; % Length of each epoch in seconds
params.EpochRest.EpochOverlap = 0.5; % Overlap between epochs (proportion)

% Epoching task data parameters (prep.segment_task)
params.EpochTask.Markers = {'instructed_toCloseEyes'}; % Event markers for epoching
params.EpochTask.TimeWindow = [1, 39]; % Time window for epoch extraction relative to marker (seconds)
