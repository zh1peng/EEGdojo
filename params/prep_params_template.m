
params=struct;

% IO parameters
% params.IO.filename='sub-NDARAD481FXF_task-RestingState_eeg.set';
% params.IO.filepath='V:\HBN_EEG_BIDS\cmi_bids_NC\sub-NDARAD481FXF\eeg';
% params.IO.outputpath='X:\HBN_EEG_preprocessing\code\test';
% params.IO = setIOPaths(params.IO);



% Channel information
params.ChanInfo.AllChan = [1:129];
params.ChanInfo.RefChan = 129;
params.ChanInfo.EOGChanIdx = [];
params.ChanInfo.EOGChanLabels = [];
params.ChanInfo.KnownBadChanIdx = []; 
params.ChanInfo.KnownBadChanLabels = [];

% Downsampling and filtering parameters
params.DownsamplingRate = 250;
params.Filter.LowCutoff = 0.1;
params.Filter.HighCutoff = 40;
params.Filter.PowerLineFrequency = 60;

% Epoching parameters
params.EpochTask.Marker = 'instructed_toCloseEyes';
params.EpochTask.EpochWindow = [1, 39];

% Bad channel detection parameters
params.BadChan.Chan2detect = setdiff(params.ChanInfo.AllChan, [params.ChanInfo.RefChan, params.ChanInfo.EOGChanIdx, params.ChanInfo.KnownBadChanIdx]);
params.BadChan.rejectchan_NormOn = 'on';
params.BadChan.rejectchan_SpecFreqRange = [0.1 40]; % Hz
params.BadChan.rejectchan_SpecThreshold = [-5 5]; % z-score
params.BadChan.rejectchan_SpecOn = 'on';
params.BadChan.rejectchan_Threshold = [-5 5]; % z-score
params.BadChan.rejectchan_KurtOn ='on';
params.BadChan.rejectchan_KurtThreshold = [-5 5]; % z-score
params.BadChan.rejectchan_ProbOn ='on';
params.BadChan.rejectchan_ProbThreshold = [-5 5]; % z-score
params.BadChan.FASTER_ChanVarOn = 'on';
params.BadChan.FASTER_ChanVarThreshold = 3; % z-score
params.BadChan.FASTER_ChanMeanCorrOn = 'on';
params.BadChan.FASTER_ChanMeanCorrThreshold = 3; % z-score
params.BadChan.FASTER_ChanHurst = 'on';
params.BadChan.FASTER_ChanHurstThreshold = 3; % z-score
params.BadChan.CleanRawData_FlatLineOn = 'on';
params.BadChan.CleanRawData_FlatLineThreshold = 5; %seconds
params.BadChan.CleanRawData_CleanDriftsHighPass = [0.25, 0.75];
params.BadChan.CleanRawData_CleanChan = 'on';
params.BadChan.CleanRawData_LineNoiseThreshold = 4;
params.BadChan.CleanRawData_CorrThreshold = 0.6;

% Re-referencing parameters
params.Refer.Chan2averageOn = 'Off';
params.Reref.Chan2average = setdiff(params.ChanInfo.AllChan, [params.ChanInfo.RefChan, params.ChanInfo.EOGChanIdx, params.ChanInfo.KnownBadChanIdx]);

% ICA parameters
params.BadIC.filterICAon = 'on';
params.BadIC.filterICAlocutoff = 1;
params.BadIC.ICAtype = 'runica';
params.BadIC.ICLabelOn = 'on';
params.BadIC.FASTEROn = 'on';
params.BadIC.ICLabelThreshold = [NaN NaN; 0.8 1; 0.8 1; 0.8 1; 0.8 1; 0.8 1; NaN NaN];
params.BadIC.Nrun = 1;

% Bad point detection parameters (clean_rawdata)
params.BadPoint.on = "on";
params.BadPoint.BurstCriterion = 20;
params.BadPoint.WindowCriterion = 0.25;
params.BadPoint.WindowCriterionTolerances = [-Inf 7];

% Epoching resting state data parameters
params.EpochRest.EpochLength = 2;
params.EpochRest.EpochOverlap = 0.5;


