classdef Params
    %PARAMS A class to hold and manage preprocessing parameters for EEGdojo.
    %   This class defines the structure for all parameters used in the
    %   EEGdojo preprocessing pipeline, ensuring consistency and easy
    %   management of settings across different processing steps.

    properties
        % General I/O and Logging
        outputpath          char
        filename            char
        filepath            char
        LogPath             char
        LogFile             char

        % Channel Information
        AllChan             double
        RefChan             double
        EOGChanIdx          double
        EOGChanLabels       cell
        KnownBadChanIdx     double
        KnownBadChanLabels  cell
        Chan2remove         {mustBeA(Chan2remove, ["double", "cell"])} % Added

        % Downsampling and Filtering
        DownsamplingRate    double
        LowCutoff           double % High-pass filter cutoff
        HighCutoff          double % Low-pass filter cutoff

        % Powerline Noise Removal
        Powerline           struct

        % Bad Channel Detection
        BadChan             struct

        % Re-referencing
        Reref               struct

        % ICA and Bad IC Detection
        BadIC               struct

        % Bad Point Detection (clean_rawdata)
        BadPoint            struct

        % Epoching Resting State Data
        EpochRest           struct

        % Epoching Task Data
        EpochTask           struct
    end

    methods
        function obj = Params(params)
            %PARAMS Construct an instance of the Params class.
            %   obj = Params(params) initializes the Params object with
            %   default values, which can be overridden by the 'params'
            %   struct provided as input.

            P = inputParser;
            P.KeepUnmatched = true; % Allow unmatched parameters to be passed through

            % General I/O and Logging
            addParameter(P, 'outputpath', pwd, @ischar);
            addParameter(P, 'filename', '', @ischar);
            addParameter(P, 'filepath', '', @ischar);
            addParameter(P, 'LogPath', pwd, @ischar);
            addParameter(P, 'LogFile', '', @ischar);

            % Channel Information
            addParameter(P, 'AllChan', [1:64], @isnumeric);
            addParameter(P, 'RefChan', 129, @isnumeric);
            addParameter(P, 'EOGChanIdx', [], @isnumeric);
            addParameter(P, 'EOGChanLabels', {}, @iscellstr);
            addParameter(P, 'KnownBadChanIdx', [], @isnumeric);
            addParameter(P, 'KnownBadChanLabels', {}, @iscellstr);
            addParameter(P, 'Chan2remove', {}, @(x) isnumeric(x) || iscellstr(x)); % Added

            % Downsampling and Filtering
            addParameter(P, 'DownsamplingRate', 250, @isnumeric);
            addParameter(P, 'LowCutoff', 0.1, @isnumeric);
            addParameter(P, 'HighCutoff', 40, @isnumeric);

            % Powerline Noise Removal (prep.remove_powerline)
            defaultPowerline = struct(...
                'Method', 'cleanline', ...
                'Freq', 60, ...
                'BW', 2, ...
                'NHarm', 3 ...
            );
            addParameter(P, 'Powerline', defaultPowerline, @isstruct);

            % Bad Channel Detection (prep.remove_bad_channels)
            defaultBadChan = struct(...
                'IdxDetect', [], ... % Will be set based on ChanInfo
                'Action', 'remove', ...
                'KnownBadIdx', [], ... % Will be set based on ChanInfo
                'Kurtosis', false, ... % Changed to false
                'Kurt_Threshold', 5, ...
                'Probability', false, ... % Changed to false
                'Prob_Threshold', 5, ...
                'Spectrum', false, ... % Changed to false
                'Spec_Threshold', 5, ...
                'Spec_FreqRange', [0.1 40], ...
                'NormOn', 'on', ...
                'FASTER_MeanCorr', false, ... % Changed to false
                'FASTER_Threshold', 0.4, ...
                'FASTER_Variance', false, ... % Changed to false
                'FASTER_VarThreshold', 3, ...
                'FASTER_Hurst', false, ... % Changed to false
                'FASTER_HurstThreshold', 3, ...
                'CleanRaw_Flatline', true, ...
                'Flatline_Sec', 5, ...
                'CleanDrift_Band', [0.25 0.75], ...
                'CleanRaw_Noise', true, ...
                'CleanChan_Corr', 0.8, ...
                'CleanChan_Line', 4, ...
                'CleanChan_MaxBad', 0.5, ...
                'CleanChan_NSamp', 50 ...
            );
            addParameter(P, 'BadChan', defaultBadChan, @isstruct);

            % Re-referencing (prep.reref)
            defaultReref = struct(...
                'excludeLabels', [] ... % Will be set based on ChanInfo
            );
            addParameter(P, 'Reref', defaultReref, @isstruct);

            % ICA and Bad IC Detection (prep.remove_bad_ICs)
            defaultBadIC = struct(...
                'FilterICAOn', true, ...
                'FilterICALocutoff', 1, ...
                'ICAType', 'runica', ...
                'ICLabelOn', true, ...
                'ICLabelThreshold', [0 0.1; 0.9 1], ...
                'FASTEROn', true, ...
                'DetectECG', false, ...
                'ECGCorrelationThreshold', 0.8 ...
            );
            addParameter(P, 'BadIC', defaultBadIC, @isstruct);

            % Bad Point Detection (clean_rawdata)
            defaultBadPoint = struct(...
                'on', 'on', ...
                'BurstCriterion', 20, ...
                'WindowCriterion', 0.25, ...
                'WindowCriterionTolerances', [-Inf 7] ...
            );
            addParameter(P, 'BadPoint', defaultBadPoint, @isstruct);

            % Epoching Resting State Data (prep.segment_rest)
            defaultEpochRest = struct(...
                'EpochLength', 2, ...
                'EpochOverlap', 0.5 ...
            );
            addParameter(P, 'EpochRest', defaultEpochRest, @isstruct);

            % Epoching Task Data (prep.segment_task)
            defaultEpochTask = struct(...
                'Markers', {'instructed_toCloseEyes'}, ...
                'TimeWindow', [1, 39] ...
            );
            addParameter(P, 'EpochTask', defaultEpochTask, @isstruct);

            % Parse the input parameters
            parse(P, params);

            % Assign parsed results to object properties
            obj.outputpath = P.Results.outputpath;
            obj.filename = P.Results.filename;
            obj.filepath = P.Results.filepath;
            obj.LogPath = P.Results.LogPath;
            obj.LogFile = P.Results.LogFile;

            obj.AllChan = P.Results.AllChan;
            obj.RefChan = P.Results.RefChan;
            obj.EOGChanIdx = P.Results.EOGChanIdx;
            obj.EOGChanLabels = P.Results.EOGChanLabels;
            obj.KnownBadChanIdx = P.Results.KnownBadChanIdx;
            obj.KnownBadChanLabels = P.Results.KnownBadChanLabels;
            obj.Chan2remove = P.Results.Chan2remove; % Added

            obj.DownsamplingRate = P.Results.DownsamplingRate;
            obj.LowCutoff = P.Results.LowCutoff;
            obj.HighCutoff = P.Results.HighCutoff;

            obj.Powerline = P.Results.Powerline;
            obj.BadChan = P.Results.BadChan;
            obj.Reref = P.Results.Reref;
            obj.BadIC = P.Results.BadIC;
            obj.BadPoint = P.Results.BadPoint;
            obj.EpochRest = P.Results.EpochRest;
            obj.EpochTask = P.Results.EpochTask;

            % Post-processing for dependent parameters
            if isempty(obj.BadChan.IdxDetect)
                obj.BadChan.IdxDetect = setdiff(obj.AllChan, [obj.RefChan, obj.EOGChanIdx, obj.KnownBadChanIdx]);
            end
            if isempty(obj.BadChan.KnownBadIdx)
                obj.BadChan.KnownBadIdx = obj.KnownBadChanIdx;
            end
            if isempty(obj.Reref.excludeLabels)
                obj.Reref.excludeLabels = setdiff(obj.AllChan, [obj.RefChan, obj.EOGChanIdx, obj.KnownBadChanIdx]);
            end
        end
    end
end