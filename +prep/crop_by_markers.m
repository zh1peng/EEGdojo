function [EEG, out] = crop_by_markers(EEG, varargin)
% CROP_BY_MARKERS Crops EEG data based on specified start and end markers.
%   This function extracts a segment of EEG data between two defined event
%   markers. An optional padding can be added before the start marker and
%   after the end marker. This is useful for isolating specific experimental
%   phases or continuous segments of interest.
%
% Syntax:
%   [EEG, out] = prep.crop_by_markers(EEG, 'param', value, ...)
%
% Input Arguments:
%   EEG         - EEGLAB EEG structure.
%
% Optional Parameters (Name-Value Pairs):
%   'StartMarker' - (char | string, default: '')
%                   The type/label of the event marker indicating the start
%                   of the segment.
%   'EndMarker'   - (char | string, default: '')
%                   The type/label of the event marker indicating the end
%                   of the segment.
%   'PadSec'      - (numeric, default: 0)
%                   Number of seconds to pad the segment on both sides.
%                   The padding is added before the start marker and after
%                   the end marker.
%   'LogFile'     - (char | string, default: '')
%                   Path to a log file for verbose output. If empty, output
%                   is directed to the command window.
%
% Output Arguments:
%   EEG         - Modified EEGLAB EEG structure containing the cropped data.
%   out         - Structure containing output information:
%                 out.start_sample - (numeric) The starting sample of the
%                                    cropped segment.
%                 out.end_sample   - (numeric) The ending sample of the
%                                    cropped segment.
%
% Examples:
%   % Example 1: Crop data between markers 'start_exp' and 'end_exp' with 1 second padding (without pipeline)
%   % Load an EEG dataset first, e.g., EEG = pop_loadset('eeg_data.set');
%   EEG_cropped = prep.crop_by_markers(EEG, ...
%       'StartMarker', 'start_exp', ...
%       'EndMarker', 'end_exp', ...
%       'PadSec', 1);
%   disp('Data cropped between "start_exp" and "end_exp" with 1s padding.');
%
%   % Example 2: Usage within a pipeline
%   % Assuming 'pipe' is an initialized pipeline object
%   pipe = pipe.addStep(@prep.crop_by_markers, ...
%       'StartMarker', 'ExperimentStart', ...
%       'EndMarker', 'ExperimentEnd', ...
%       'PadSec', 0.5, ...
%       'LogFile', p.LogFile); %% p.LogFile from pipeline parameters
%   % Then run the pipeline: [EEG_processed, results] = pipe.run(EEG);
%   disp('Data cropped via pipeline.');
%
% See also: pop_select, eeg_checkset

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('StartMarker', '', @(s) ischar(s) || isstring(s)); % Allow string type
    p.addParameter('EndMarker', '', @(s) ischar(s) || isstring(s));   % Allow string type
    p.addParameter('PadSec', 0, @isnumeric);
    p.addParameter('LogFile', '', @(s) ischar(s) || isstring(s));

    p.parse(EEG, varargin{:});
    R = p.Results;

    % Extract parameters for clarity and ensure char type for comparison
    StartMarker = char(R.StartMarker);
    EndMarker = char(R.EndMarker);
    PadTime = R.PadSec;
    LogFile = R.LogFile;

    out = struct(); % Initialize output structure


    if isempty(StartMarker) || isempty(EndMarker)
        error('[crop_by_markers] Both StartMarker and EndMarker must be specified.');
    end

    logPrint(LogFile, sprintf('[crop_by_markers] Cropping EEG data between markers "%s" and "%s" with %.2f seconds padding.', StartMarker, EndMarker, PadTime));

    % Ensure event types are consistent (char)
    evtype = {EEG.event.type};
    for i = 1:numel(evtype)
        if ~(ischar(evtype{i}) || isstring(evtype{i}))
            evtype{i} = num2str(evtype{i}); % Convert numeric types to char
        else
            evtype{i} = char(evtype{i});    % Ensure string types are char
        end
    end

    % Find start marker
    start_idx = find(strcmp(evtype, StartMarker), 1, 'first');
    if isempty(start_idx)
        error('[crop_by_markers] Start marker "%s" not found in EEG events.', StartMarker);
    end

    % Find end marker (must occur after start marker)
    lat_all   = [EEG.event.latency];
    end_all   = find(strcmp(evtype, EndMarker));
    end_after = end_all(lat_all(end_all) > lat_all(start_idx));
    if isempty(end_after)
        error('[crop_by_markers] End marker "%s" not found after start marker "%s".', EndMarker, StartMarker);
    end
    end_idx = end_after(end); % Use the last occurrence after the start marker


    % Calculate start and end samples for cropping
    start_samp = EEG.event(start_idx).latency;
    end_samp   = EEG.event(end_idx).latency;

    pad_samp  = round(max(0, PadTime) * EEG.srate); % Convert padding seconds to samples
    seg_start = max(1, start_samp - pad_samp);     % Ensure start is not before 1
    seg_end   = min(EEG.pnts, end_samp + pad_samp); % Ensure end is not after total points

    % Perform the cropping
    EEG = pop_select(EEG, 'point', [seg_start, seg_end]);
    EEG = eeg_checkset(EEG); % Update EEG structure after changes

    logPrint(LogFile, sprintf('[crop_by_markers] Data successfully cropped from sample %d to %d.', seg_start, seg_end));
    
    % Store output information
    out.start_sample = seg_start;
    out.end_sample = seg_end;

end