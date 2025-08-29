function [EEG, out] = interpolate_bad_channels_epoch(EEG, varargin)
% INTERPOLATE_BAD_CHANNELS_EPOCH  Identifies and interpolates bad channels at the epoch level.
%   This function uses the FASTER algorithm to find bad channels for each epoch
%   and then interpolates them using spherical splines.
%
% Syntax:
%   [EEG, out] = prep.interpolate_bad_channels_epoch(EEG, 'param', value, ...)
%
% Input Arguments:
%   EEG         - EEGLAB EEG structure (epoched data).
%
% Optional Parameters (Name-Value Pairs):
%   'LogFile'           - (char | string, default: '')
%                         File path to log the results.
%
% Output Arguments:
%   EEG         - Modified EEGLAB EEG structure with bad channels interpolated.
%   out         - Structure containing details of the detection:
%                 out.bad_chan_cell: Cell array with bad channel indices for each epoch.
%
% See also: single_epoch_channel_properties, h_epoch_interp_spl

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('LogFile', '', @(s) ischar(s) || isstring(s));
    p.parse(EEG, varargin{:});
    R = p.Results;

    out = struct(); % Initialize out struct

    logPrint(R.LogFile, 'Identifying bad channels per epoch...');

    % ----------------- Find bad channels per epoch -----------------
    bad_chan_cell = cell(1, EEG.trials);
    for epoch_i = 1:EEG.trials
        bad_chan_epoch_list = single_epoch_channel_properties(EEG, epoch_i, 1:EEG.nbchan);
        tmp_bad = find(min_z(bad_chan_epoch_list) == 1);
        bad_chan_cell{epoch_i} = tmp_bad;
        if ~isempty(tmp_bad)
            logPrint(R.LogFile, sprintf('Epoch %d - Bad Channels: %d, Details: %s', epoch_i, length(tmp_bad), mat2str(tmp_bad)));
        end
    end

    % ----------------- Interpolate bad channels -----------------
    logPrint(R.LogFile, 'Interpolating bad channels at the epoch level...');
    EEG = h_epoch_interp_spl(EEG, bad_chan_cell);
    logPrint(R.LogFile, 'Bad channels interpolated successfully.');

    % ----------------- Bookkeeping in EEG.etc -----------------
    out.bad_chan_cell = bad_chan_cell;

    if ~isfield(EEG.etc, 'EEGdojo'), EEG.etc.EEGdojo = struct(); end
    EEG.etc.EEGdojo.BadChanEpochCell = bad_chan_cell;

end