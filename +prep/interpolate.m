function [EEG, out] = interpolate(EEG, varargin)
% INTERPOLATE Interpolates missing channels in EEG data.
%
% This function identifies channels present in the original channel locations
% (`EEG.urchanlocs`) but missing from the current EEG dataset (`EEG.chanlocs`),
% and then interpolates their data using spherical interpolation. This is
% typically used to restore data for channels that were removed (e.g., bad channels).
%
% Inputs:
%   EEG         - EEGLAB EEG structure. Must contain `EEG.urchanlocs` with
%                 original channel locations.
%   varargin    - Optional parameters:
%     'LogFile' - (char) Path to the log file for recording processing
%                 information. Default is ''.
%
% Outputs:
%   EEG         - Modified EEGLAB EEG structure with interpolated channels.
%   out         - Structure containing output information:
%     .interpolated_channels - (cell array of strings) Labels of channels that were interpolated.
%
% Examples:
%   % 1. Interpolate missing channels:
%   EEG = interpolate(EEG);
%
%   % 2. Usage within a pipeline:
%   %    (Assuming 'p' is a parameter structure containing 'p.logFile')
%   pipe = pipe.addStep(@prep.interpolate, ...
%       'LogFile', p.logFile);
%
% See also: pop_interp, eeg_checkset

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('LogFile', '', @ischar);

    p.parse(EEG, varargin{:});
    R = p.Results;

    out = struct(); % Initialize output structure
    out.interpolated_channels = {};

    try
        logPrint(R.LogFile, '[interpolate] Starting channel interpolation.');

        if ~isfield(EEG, 'urchanlocs') || isempty(EEG.urchanlocs)
            error('[interpolate] Skipping interpolation: No original channel locations (EEG.urchanlocs) found.');
        end

        original_chans = {EEG.urchanlocs.labels};
        current_chans = {EEG.chanlocs.labels};
        chans_to_interp = setdiff(original_chans, current_chans);

        if isempty(chans_to_interp)
           logPrint(R.LogFile, '[interpolate] Skipping interpolation: No channels to interpolate found.');
           return;
        end

        logPrint(R.LogFile, sprintf('[interpolate] Interpolating %d channels: %s', numel(chans_to_interp), strjoin(chans_to_interp, ', ')));
        EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');
        EEG = eeg_checkset(EEG); % Update EEG structure after changes
        logPrint(R.LogFile, '[interpolate] Interpolation complete.');
        out.interpolated_channels = chans_to_interp;

    catch ME
        error('[interpolate] Interpolation failed: %s', ME.message);
    end
end
