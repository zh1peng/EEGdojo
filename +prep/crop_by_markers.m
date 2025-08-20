function EEG = crop_by_markers(EEG, varargin)
%CROP_BY_MARKERS Crop EEG data based on start and end markers.
%   Usage: EEG = crop_by_markers(EEG, 'StartMarker','1', 'EndMarker','2', 'PadSec', pad_seconds);

    p = inputParser;
    addRequired(p, 'EEG', @isstruct);
    addRequired(p, 'StartMarker', @ischar);
    addRequired(p, 'EndMarker', @ischar);
    addParameter(p, 'PadSec', 0, @isnumeric);
    addParameter(p, 'logFile', '', @ischar);
    parse(p, EEG, varargin{:});

    
    R = p.Results;
    EEG = R.EEG;
    StartMarker = R.StartMarker;
    EndMarker = R.EndMarker;
    PadTime = R.PadSec;

    logPrint(R.logFile, sprintf('--- Crop EEG data between markers "%s" and "%s" ---', StartMarker, EndMarker));
    evtype = {EEG.event.type};
    for i = 1:numel(evtype)
        if ~(ischar(evtype{i}) || isstring(evtype{i}))
            evtype{i} = num2str(evtype{i});
        else
            evtype{i} = char(evtype{i});
        end
    end

    start_idx = find(strcmp(evtype, StartMarker), 1, 'first');
    if isempty(start_idx)
        error('Start marker "%s" not found.', StartMarker);
    end

    lat_all   = [EEG.event.latency];
    end_all   = find(strcmp(evtype, EndMarker));
    end_after = end_all(lat_all(end_all) > lat_all(start_idx));
    if isempty(end_after)
        error('End marker "%s" not found.', EndMarker);
    end
    end_idx = end_after(end);

    try 
        start_samp = EEG.event(start_idx).latency;
        end_samp   = EEG.event(end_idx).latency;

        pad_samp  = round(max(0, PadTime) * EEG.srate);
        seg_start = max(1, start_samp - pad_samp);
        seg_end   = min(EEG.pnts, end_samp + pad_samp);

        EEG = pop_select(EEG, 'point', [seg_start, seg_end]);
        EEG = eeg_checkset(EEG);
        logPrint(R.logFile, sprintf('Segmenting data from %s to %s.', EEG.event(start_idx).label, EEG.event(end_idx).label));
    catch ME 
        error('Segmentation failed: %s.',ME.message);
    end
end
