classdef Params
    %PARAMS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        outputpath
        filename
        filepath
        StartMarker
        EndMarker
        PadTime
        eogchan
        ICLabel
        chan2remove
        DownsamplingRate
        HPfreq
        LPfreq
        PowerLineFrequency
        FlatLineCriterion
        ChannelCriterion
        LineNoiseCriterion
        BadChan
    end

    methods
        function obj = Params(params)
            %PARAMS Construct an instance of this class
            %   Detailed explanation goes here
            P = inputParser;
            P.KeepUnmatched = true;
            
            addParameter(P, 'outputpath', pwd, @ischar);
            addParameter(P, 'filename', '', @ischar);
            addParameter(P, 'filepath', '', @ischar);
            addParameter(P, 'StartMarker', '', @ischar);
            addParameter(P, 'EndMarker', '', @ischar);
            addParameter(P, 'PadTime', 0, @isnumeric);
            addParameter(P, 'eogchan', {}, @(x) iscell(x) || ischar(x));
            addParameter(P, 'ICLabel', [0.6 1; 0.6 1; 0.6 1], @isnumeric);
            addParameter(P, 'chan2remove', {}, @(x) iscell(x) || ischar(x));
            addParameter(P, 'DownsamplingRate', 250, @isnumeric);
            addParameter(P, 'HPfreq', 0.5, @isnumeric);
            addParameter(P, 'LPfreq', 40, @isnumeric);
            addParameter(P, 'PowerLineFrequency', 60, @isnumeric);
            addParameter(P, 'FlatLineCriterion', 5, @isnumeric);
            addParameter(P, 'ChannelCriterion', 0.8, @isnumeric);
            addParameter(P, 'LineNoiseCriterion', 4, @isnumeric);

            defaultBadChan = struct(...
                'rejectchan_KurtOn', 'off', ...
                'rejectchan_KurtThreshold', 5, ...
                'rejectchan_NormOn', 'on', ...
                'rejectchan_SpecOn', 'off', ...
                'rejectchan_SpecThreshold', 5, ...
                'rejectchan_SpecFreqRange', [1 50], ...
                'rejectchan_ProbOn', 'off', ...
                'rejectchan_ProbThreshold', 5, ...
                'RemoveRefChan', false, ...
                'RefChan', [] ...
            );
            addParameter(P, 'BadChan', defaultBadChan, @isstruct);
            
            parse(P, params);
            
            obj.outputpath = P.Results.outputpath;
            obj.filename = P.Results.filename;
            obj.filepath = P.Results.filepath;
            obj.StartMarker = P.Results.StartMarker;
            obj.EndMarker = P.Results.EndMarker;
            obj.PadTime = P.Results.PadTime;
            obj.eogchan = P.Results.eogchan;
            obj.ICLabel = P.Results.ICLabel;
            obj.chan2remove = P.Results.chan2remove;
            obj.DownsamplingRate = P.Results.DownsamplingRate;
            obj.HPfreq = P.Results.HPfreq;
            obj.LPfreq = P.Results.LPfreq;
            obj.PowerLineFrequency = P.Results.PowerLineFrequency;
            obj.FlatLineCriterion = P.Results.FlatLineCriterion;
            obj.ChannelCriterion = P.Results.ChannelCriterion;
            obj.LineNoiseCriterion = P.Results.LineNoiseCriterion;
            obj.BadChan = P.Results.BadChan;
        end
    end
end
