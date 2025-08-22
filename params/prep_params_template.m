Params=struct;

%% General I/O and Logging parameters
Params.IO.filename=''; % Example: 'sub-NDARAD481FXF_task-RestingState_eeg.set';
Params.IO.filepath=''; % Example: 'V:\HBN_EEG_BIDS\cmi_bids_NC\sub-NDARAD481FXF\eeg';
Params.IO.outputpath=''; % Example: 'X:\HBN_EEG_preprocessing\code\test';

%% Channel information we only use label as using idx is very risky when there are removel of channels

Params.ChanInfo.RefChanLabels = {}; % Reference channel index (if applicable)
Params.ChanInfo.EOGChanLabels = {'VEOU','VEOL','HEOR','HEOL'}; % Labels of EOG channels (e.g., {'VEOG', 'HEOG'})
Params.ChanInfo.KnownBadChanLabels = {}; % Labels of channels known to be bad
Params.ChanInfo.Chan2remove = {}; % Channels to remove (e.g., {'ref''} 
Params.ChanInfo.ECGChanLabel = 'ECG';
Params.ChanInfo.OtherChanLabel = {};

%% continuse data
Params.Crop.StartMarker = '1'; 
Params.Crop.EndMarker   = '2';
Params.Crop.PadTime     = 0.5;  

%% Downsampling and filtering parameters
Params.DownsamplingRate = 250; % Target downsampling rate in Hz
Params.Filter.LowCutoff = 0.1; % High-pass filter cutoff in Hz
Params.Filter.HighCutoff = 40; % Low-pass filter cutoff in Hz

% Powerline noise removal parameters (prep.remove_powerline)
Params.Powerline.Method = 'cleanline'; % 'cleanline' or 'notch'
Params.Powerline.Freq = 50; % Fundamental powerline frequency (e.g., 50 or 60 Hz)
Params.Powerline.BW = 2; % Half-bandwidth for FIR notch filter (only for 'notch' method)
Params.Powerline.NHarm = 3; % Number of harmonics to target

% Bad channel detection parameters (prep.remove_bad_channels)
Params.BadChan.ExcludeLabel = [Params.ChanInfo.EOGChanLabels, Params.ChanInfo.RefChanLabels, Params.ChanInfo.OtherChanLabel, Params.ChanInfo.ECGChanLabel]; % Channels to consider for detection
Params.BadChan.Action = 'remove'; % 'remove' or 'flag'
Params.BadChan.KnownBadLabel = Params.ChanInfo.KnownBadChanLabels; % Indices of channels already known to be bad

% Classic EEGLAB detectors (pop_rejchan)
Params.BadChan.Kurtosis = false; % Enable kurtosis-based detection
Params.BadChan.Kurt_Threshold = 5; % Kurtosis threshold (z-score)
Params.BadChan.Probability = false; % Enable probability-based detection
Params.BadChan.Prob_Threshold = 5; % Probability threshold (z-score)
Params.BadChan.Spectrum = false; % Enable spectrum-based detection
Params.BadChan.Spec_Threshold = 5; % Spectrum threshold (z-score)
Params.BadChan.Spec_FreqRange = [0.1 40]; % Frequency range for spectrum analysis (Hz)
Params.BadChan.NormOn = 'on'; % Normalize measures ('on' or 'off')

% FASTER detectors (FASTER_rejchan)
Params.BadChan.FASTER_MeanCorr = false; % Enable FASTER mean correlation
Params.BadChan.FASTER_Threshold = 0.4; % FASTER mean correlation threshold
Params.BadChan.FASTER_Variance = false; % Enable FASTER variance
Params.BadChan.FASTER_VarThreshold = 3; % FASTER variance threshold
Params.BadChan.FASTER_Hurst = false; % Enable FASTER Hurst exponent
Params.BadChan.FASTER_HurstThreshold = 3; % FASTER Hurst exponent threshold

% CleanRaw detectors (cleanraw_rejchan)
Params.BadChan.CleanRaw_Flatline = true; % Enable CleanRaw flatline detection
Params.BadChan.Flatline_Sec = 5; % Minimum flatline duration in seconds
Params.BadChan.CleanDrift_Band = [0.25 0.75]; % Highpass filter band for CleanRaw methods
Params.BadChan.CleanRaw_Noise = true; % Enable CleanRaw noise detection (channel correlation)
Params.BadChan.CleanChan_Corr = 0.8; % Correlation threshold for CleanRaw noise
Params.BadChan.CleanChan_Line = 4; % Line noise threshold for CleanRaw noise
Params.BadChan.CleanChan_MaxBad = 0.5; % Maximum proportion of bad time points for CleanRaw noise
Params.BadChan.CleanChan_NSamp = 50; % Number of samples for CleanRaw noise

%% Re-referencing parameters (prep.reref)
Params.Reref.excludeLabels = {}; % Channels to exclude from average reference

%% ICA parameters (prep.remove_bad_ICs)
Params.BadIC.FilterICAOn = true; % Apply high-pass filter before ICA
Params.BadIC.FilterICALocutoff = 1; % High-pass filter cutoff for ICA (Hz)
Params.BadIC.ICAType = 'runica'; % Type of ICA algorithm
Params.BadIC.ICLabelOn = true; % Enable ICLabel classification
Params.BadIC.ICLabelThreshold = [NaN NaN; 0.7 1; 0.7 1; 0.7 1; 0.7 1; 0.7 1; NaN NaN]; % ICLabel thresholds
Params.BadIC.FASTEROn = true; % Enable FASTER component property analysis
Params.BadIC.FASTERExludeTreshold = 0.7; % Exclude components with proablity of being brain signal by IClabel > threshold
Params.BadIC.DetectECG = false; % Enable ECG correlation-based detection
Params.BadIC.ECGStruct = []; % New: Placeholder for external ECG data structure
Params.BadIC.ECGCorrelationThreshold = 0.8; % ECG correlation threshold



% % Epoching resting state data parameters (prep.segment_rest)
% Params.EpochRest.EpochLength = 2; % Length of each epoch in seconds
% Params.EpochRest.EpochOverlap = 0.5; % Overlap between epochs (proportion)

% % Epoching task data parameters (prep.segment_task)
% Params.EpochTask.Markers = {'instructed_toCloseEyes'}; % Event markers for epoching
% Params.EpochTask.TimeWindow = [1, 39]; % Time window for epoch extraction relative to marker (seconds)



