function [EEG, out] = remove_bad_channels(EEG, varargin)
% REMOVE_BAD_CHANNELS  Detects and removes/flags bad EEG channels using various methods.
%   This function provides a modular approach to identify bad channels based on
%   EEGLAB's `pop_rejchan` (kurtosis, probability, spectrum), FASTER (mean
%   correlation, variance, Hurst exponent), and CleanRaw (flatline, noise)
%   algorithms. It can either remove the identified bad channels from the EEG
%   dataset or flag them by adding a mask to `EEG.etc.clean_channel_mask`.
%
% Syntax:
%   [EEG, out] = prep.remove_bad_channels(EEG, 'param', value, ...)
%
% Input Arguments:
%   EEG         - EEGLAB EEG structure.
%
% Optional Parameters (Name-Value Pairs):
%   'IdxDetect'         - (numeric array | logical array, default: all channels)
%                         Indices of channels to consider for detection. If empty,
%                         all channels are considered.
%   'Action'            - (char | string, 'remove' | 'flag', default: 'remove')
%                         'remove': Removes bad channels from the dataset.
%                         'flag': Adds a mask to EEG.etc.clean_channel_mask
%                                 indicating good (true) and bad (false) channels.
%   'LogPath'           - (char | string, default: pwd)
%                         Path to save log plots and reports.
%   'LogFile'           - (char | string, default: '')
%                         Base name for the log report file.
%   'KnownBadIdx'       - (numeric array, default: [])
%                         Indices of channels already known to be bad. These will
%                         be included in the final list of bad channels.
%
%   %% Classic EEGLAB detectors (pop_rejchan)
%   'Kurtosis'          - (logical, default: false)
%                         Enable kurtosis-based channel rejection.
%   'Kurt_Threshold'    - (numeric, default: 5)
%                         Threshold for kurtosis.
%   'Probability'       - (logical, default: false)
%                         Enable probability-based channel rejection.
%   'Prob_Threshold'    - (numeric, default: 5)
%                         Threshold for probability.
%   'Spectrum'          - (logical, default: false)
%                         Enable spectrum-based channel rejection.
%   'Spec_Threshold'    - (numeric, default: 5)
%                         Threshold for spectrum.
%   'Spec_FreqRange'    - (numeric array [min_freq max_freq], default: [1 50])
%                         Frequency range for spectrum analysis.
%   'NormOn'            - (char | string, 'on' | 'off', default: 'on')
%                         Normalize measures ('on' or 'off') for pop_rejchan.
%
%   %% FASTER detectors (FASTER_rejchan)
%   'FASTER_MeanCorr'   - (logical, default: false)
%                         Enable FASTER mean correlation rejection.
%   'FASTER_Threshold'  - (numeric, default: 0.4)
%                         Threshold for FASTER mean correlation.
%   'FASTER_RefChan'    - (numeric, default: [])
%                         Reference channel index for FASTER methods.
%   'FASTER_Bandpass'   - (numeric array [low_freq high_freq], default: [])
%                         Bandpass filter for FASTER methods.
%   'FASTER_Variance'   - (logical, default: false)
%                         Enable FASTER variance rejection.
%   'FASTER_VarThreshold'- (numeric, default: 3)
%                         Threshold for FASTER variance.
%   'FASTER_Hurst'      - (logical, default: false)
%                         Enable FASTER Hurst exponent rejection.
%   'FASTER_HurstThreshold'- (numeric, default: 3)
%                         Threshold for FASTER Hurst exponent.
%
%   %% CleanRaw detectors (cleanraw_rejchan)
%   'CleanRaw_Flatline' - (logical, default: false)
%                         Enable CleanRaw flatline rejection.
%   'Flatline_Sec'      - (numeric, default: 5)
%                         Minimum flatline duration in seconds.
%   'CleanDrift_Band'   - (numeric array [low_freq high_freq], default: [0.25 0.75])
%                         Highpass filter band for CleanRaw methods.
%   'CleanRaw_Noise'    - (logical, default: false)
%                         Enable CleanRaw noise rejection (channel correlation).
%   'CleanChan_Corr'    - (numeric, default: 0.8)
%                         Correlation threshold for CleanRaw noise.
%   'CleanChan_Line'    - (numeric, default: 4)
%                         Line noise threshold for CleanRaw noise.
%   'CleanChan_MaxBad'  - (numeric, default: 0.5)
%                         Maximum proportion of bad time points for CleanRaw noise.
%   'CleanChan_NSamp'   - (numeric, default: 50)
%                         Number of samples for CleanRaw noise.
%
% Output Arguments:
%   EEG         - Modified EEGLAB EEG structure with bad channels removed or flagged.
%   out         - Structure containing details of the detection:
%                 out.Bad: Structure with indices of bad channels per detector.
%                 out.Bad.all: Combined unique indices of all bad channels.
%                 out.summary: Summary of bad channels per detector and total.
%                 out.detectors_used: Cell array of detector names used.
%                 out.IdxDetect: Indices of channels considered for detection.
%
% Examples:
%   % Example 1: Remove bad channels using Kurtosis and Probability (without pipeline)
%   % Load an EEG dataset first, e.g., EEG = pop_loadset('eeg_data.set');
%   [EEG_cleaned, bad_chan_info] = prep.remove_bad_channels(EEG, ...
%       'Kurtosis', true, 'Kurt_Threshold', 5, ...
%       'Probability', true, 'Prob_Threshold', 5, ...
%       'LogPath', 'C:\temp\eeg_logs', 'LogFile', 'bad_channels_report');
%   disp('Bad channels removed:');
%   disp(bad_chan_info.Bad.all);
%
%   % Example 2: Flag bad channels using FASTER Mean Correlation (with pipeline)
%   % Assuming 'pipe' is an initialized pipeline object
%   pipe = pipe.addStep(@prep.remove_bad_channels, ...
%       'Action', 'flag', ...
%       'FASTER_MeanCorr', true, 'FASTER_Threshold', 0.3, ...
%       'LogPath', 'C:\temp\eeg_logs', 'LogFile', 'flagged_channels_report');
%   % Then run the pipeline: [EEG_processed, results] = pipe.run(EEG);
%   % The flagged channels will be in EEG_processed.etc.clean_channel_mask
%
% See also: pop_rejchan, FASTER_rejchan, cleanraw_rejchan, logplot_badchannels, logreport_badchannels

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);

    % Scope, I/O, and action parameters
    p.addParameter('IdxDetect',         [], @(x) isnumeric(x) || islogical(x));
    p.addParameter('Action',            'remove', @(s) any(strcmpi(s,{'remove','flag'})));
    p.addParameter('LogPath',           '', @(s) ischar(s) || isstring(s));
    p.addParameter('LogFile',           '', @(s) ischar(s) || isstring(s));
    p.addParameter('KnownBadIdx',       [], @(x) isempty(x) || isnumeric(x));

    % Classic EEGLAB detectors
    p.addParameter('Kurtosis',          false, @islogical);
    p.addParameter('Kurt_Threshold',    5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    p.addParameter('Probability',       false, @islogical);
    p.addParameter('Prob_Threshold',    5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    p.addParameter('Spectrum',          false, @islogical);
    p.addParameter('Spec_Threshold',    5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    p.addParameter('Spec_FreqRange',    [1 50], @(x)isnumeric(x)&&numel(x)==2&&all(x>=0));
    p.addParameter('NormOn',            'on', @(x) any(strcmpi(x,{'on','off'})));

    % FASTER detectors
    p.addParameter('FASTER_MeanCorr',     false, @islogical);
    p.addParameter('FASTER_Threshold',    0.4, @isnumeric);
    p.addParameter('FASTER_RefChan',      [], @(x) isempty(x) || (isscalar(x)&&isnumeric(x)));
    p.addParameter('FASTER_Bandpass',     [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
    p.addParameter('FASTER_Variance',     false, @islogical);
    p.addParameter('FASTER_VarThreshold', 3, @isnumeric);
    p.addParameter('FASTER_Hurst',        false, @islogical);
    p.addParameter('FASTER_HurstThreshold', 3, @isnumeric);

    % CleanRaw detectors
    p.addParameter('CleanRaw_Flatline', false, @islogical);
    p.addParameter('Flatline_Sec',      5, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
    p.addParameter('CleanDrift_Band',   [0.25 0.75], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
    p.addParameter('CleanRaw_Noise',    false, @islogical);
    p.addParameter('CleanChan_Corr',    0.8, @isnumeric);
    p.addParameter('CleanChan_Line',    4, @isnumeric);
    p.addParameter('CleanChan_MaxBad',  0.5, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
    p.addParameter('CleanChan_NSamp',   50, @isnumeric);

    p.parse(EEG, varargin{:});
    R = p.Results;

    out = struct(); % Initialize out struct


    % ----------------- Initial Setup -----------------
    if isempty(R.IdxDetect), R.IdxDetect = 1:EEG.nbchan; end
    if islogical(R.IdxDetect), R.IdxDetect = find(R.IdxDetect); end
    if ~exist(R.LogPath,'dir')&~isempty(R.LogPath), mkdir(R.LogPath); end

    % Preserve original channel locations
    if ~isfield(EEG,'urchanlocs') || isempty(EEG.urchanlocs), [EEG.urchanlocs] = deal(EEG.chanlocs); end
    if ~isfield(EEG,'chaninfo') || ~isfield(EEG.chaninfo,'nodatchans'), EEG.chaninfo.nodatchans = []; end

    % ----------------- Run Detectors -----------------
    Bad  = struct();
    used = {};

    if R.Kurtosis
        logPrint(R.LogFile, '[remove_bad_channels] Running Kurtosis detector...');
        [~, idxRel] = pop_rejchan(EEG, 'elec', R.IdxDetect, 'threshold', R.Kurt_Threshold, 'norm', lower(R.NormOn), 'measure', 'kurt');
        Bad.Kurt = R.IdxDetect(idxRel);
        used{end+1} = 'Kurt';
        logplot_badchannels(EEG, Bad.Kurt, R.LogPath, 'Kurt');
    else
        Bad.Kurt = [];
    end

    if R.Spectrum
        logPrint(R.LogFile, '[remove_bad_channels] Running Spectrum detector...');
        [~, idxRel] = pop_rejchan(EEG, 'elec', R.IdxDetect, 'threshold', R.Spec_Threshold, 'norm', lower(R.NormOn), 'measure', 'spec', 'freqrange', R.Spec_FreqRange);
        Bad.Spec = R.IdxDetect(idxRel);
        used{end+1} = 'Spec';
        logplot_badchannels(EEG, Bad.Spec, R.LogPath, 'Spec');
    else
        Bad.Spec = [];
    end

    if R.Probability
        logPrint(R.LogFile, '[remove_bad_channels] Running Probability detector...');
        [~, idxRel] = pop_rejchan(EEG, 'elec', R.IdxDetect, 'threshold', R.Prob_Threshold, 'norm', lower(R.NormOn), 'measure', 'prob');
        Bad.Prob = R.IdxDetect(idxRel);
        used{end+1} = 'Prob';
        logplot_badchannels(EEG, Bad.Prob, R.LogPath, 'Prob');
    else
        Bad.Prob = [];
    end

    if R.FASTER_MeanCorr
        logPrint(R.LogFile, '[remove_bad_channels] Running FASTER Mean Correlation detector...');
        args = {'elec', R.IdxDetect, 'measure','meanCorr', 'threshold', R.FASTER_Threshold};
        if ~isempty(R.FASTER_RefChan), args = [args, {'refchan', R.FASTER_RefChan}]; end
        if ~isempty(R.FASTER_Bandpass), args = [args, {'bandpass', R.FASTER_Bandpass}]; end
        [~, idx] = FASTER_rejchan(EEG, args{:});
        if ~isempty(R.FASTER_RefChan), idx = setdiff(idx, R.FASTER_RefChan); end
        Bad.MeanCorr = idx(:)';
        used{end+1} = 'FASTER_MeanCorr';
        logplot_badchannels(EEG, Bad.MeanCorr, R.LogPath, 'MeanCorr');
    else
        Bad.MeanCorr = [];
    end

    if R.FASTER_Variance
        logPrint(R.LogFile, '[remove_bad_channels] Running FASTER Variance detector...');
        args = {'elec', R.IdxDetect, 'measure','variance', 'threshold', R.FASTER_VarThreshold};
        if ~isempty(R.FASTER_RefChan), args = [args, {'refchan', R.FASTER_RefChan}]; end
        if ~isempty(R.FASTER_Bandpass), args = [args, {'bandpass', R.FASTER_Bandpass}]; end
        [~, idx] = FASTER_rejchan(EEG, args{:});
        if ~isempty(R.FASTER_RefChan), idx = setdiff(idx, R.FASTER_RefChan); end
        Bad.Variance = idx(:)';
        used{end+1} = 'FASTER_Variance';
        logplot_badchannels(EEG, Bad.Variance, R.LogPath, 'Variance');
    else
        Bad.Variance = [];
    end

    if R.FASTER_Hurst
        logPrint(R.LogFile, '[remove_bad_channels] Running FASTER Hurst detector...');
        args = {'elec', R.IdxDetect, 'measure','hurst', 'threshold', R.FASTER_HurstThreshold};
        if ~isempty(R.FASTER_RefChan), args = [args, {'refchan', R.FASTER_RefChan}]; end
        if ~isempty(R.FASTER_Bandpass), args = [args, {'bandpass', R.FASTER_Bandpass}]; end
        [~, idx] = FASTER_rejchan(EEG, args{:});
        if ~isempty(R.FASTER_RefChan), idx = setdiff(idx, R.FASTER_RefChan); end
        Bad.Hurst = idx(:)';
        used{end+1} = 'FASTER_Hurst';
        logplot_badchannels(EEG, Bad.Hurst, R.LogPath, 'Hurst');
    else
            Bad.Hurst = [];
    end

    if R.CleanRaw_Flatline
        logPrint(R.LogFile, '[remove_bad_channels] Running CleanRaw Flatline detector...');
        idx = cleanraw_rejchan(EEG, 'elec', R.IdxDetect, 'measure','flatline', 'threshold', R.Flatline_Sec, 'highpass', R.CleanDrift_Band);
        Bad.Flatline = idx(:)';
        used{end+1} = 'Flatline';
        logplot_badchannels(EEG, Bad.Flatline, R.LogPath, 'Flatline');
    else
        Bad.Flatline = [];
    end

    if R.CleanRaw_Noise
        logPrint(R.LogFile, '[remove_bad_channels] Running CleanRaw Noise detector...');
        idx = cleanraw_rejchan(EEG, 'elec', R.IdxDetect, 'measure','CleanChan', 'chancorr_crit', R.CleanChan_Corr, 'line_crit', R.CleanChan_Line, 'max_broken_time', R.CleanChan_MaxBad, 'min_corr_samples', R.CleanChan_NSamp, 'highpass', R.CleanDrift_Band);
        Bad.CleanChan = idx(:)';
        used{end+1} = 'CleanChan';
        logplot_badchannels(EEG, Bad.CleanChan, R.LogPath, 'CleanChan');
    else
        Bad.CleanChan = [];
    end

    Bad.Known = unique(R.KnownBadIdx(:)','stable');

    % ----------------- Combine & Summarize -----------------
    fields = {'Kurt','Spec','Prob','MeanCorr','Variance','Hurst','Flatline','CleanChan','Known'};
    allBad = [];
    summary = struct();
    for k = 1:numel(fields)
        f = fields{k};
        if ~isfield(Bad,f) || isempty(Bad.(f)), Bad.(f) = []; end
        summary.(f) = numel(Bad.(f));
        allBad = [allBad, Bad.(f)]; %#ok<AGROW>
    end
    Bad.all = unique(allBad, 'stable');
    summary.Total = numel(Bad.all);

    logPrint(R.LogFile, sprintf('[remove_bad_channels] Total bad channels identified: %d', summary.Total));

    % ----------------- Log Report -----------------
    logreport_badchannels(Bad, R.LogFile);

    % ----------------- Action: Remove or Flag -----------------
    switch lower(R.Action)
        case 'remove'
            if ~isempty(Bad.all)
                logPrint(R.LogFile,sprintf( '[remove_bad_channels] Removing %d bad channels...', summary.Total));
                EEG = pop_select(EEG, 'rmchannel', Bad.all);
                EEG = eeg_checkset(EEG);
                logPrint(R.LogFile, '[remove_bad_channels] Bad channels removed successfully.');
            else
                logPrint(R.LogFile, '[remove_bad_channels] No bad channels to remove.');
            end
        case 'flag'
            logPrint(R.LogFile, sprintf('[remove_bad_channels] Flagging %d bad channels...', summary.Total));
            mask = true(1, EEG.nbchan);
            mask(Bad.all(Bad.all>=1 & Bad.all<=EEG.nbchan)) = false;
            EEG.etc.clean_channel_mask = mask(:)';
            logPrint(R.LogFile, '[remove_bad_channels] Bad channels flagged successfully.');
    end

    % ----------------- Bookkeeping in EEG.etc -----------------
    out.Bad = Bad;
    out.summary = summary;
    out.detectors_used = used;
    out.IdxDetect = R.IdxDetect;

    if ~isfield(EEG.etc, 'EEGdojo'), EEG.etc.EEGdojo = struct(); end
    EEG.etc.EEGdojo.BadChanIdx = Bad.all;
    EEG.etc.EEGdojo.BadChanLabel = idx2chans(EEG, Bad.all);
    EEG.etc.EEGdojo.BadChanSummary = summary;
    EEG.etc.EEGdojo.BadDetectorsUsed = used;
end