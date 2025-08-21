function [EEG, out] = load_set(~, varargin)
% LOAD_SET  Load an EEGLAB .set dataset
%
% Usage:
%   [EEG, out] = eegdojo.io.load_set([], 'filename','x.set','filepath','./data')
%
    p = inputParser;
    p.addParameter('filename','',@ischar);
    p.addParameter('filepath','',@ischar); % Changed from addRequired
    p.addParameter('logFile', '', @ischar); % Changed from addRequired
    p.parse(varargin{:});

    R = p.Results;
    logPrint(R.logFile, sprintf('[load_set] --- Loading dataset: %s/%s ---', R.filepath, R.filename)); % Added [load_set] prefix

    EEG = pop_loadset('filename', R.filename, 'filepath', R.filepath);
    EEG = eeg_checkset(EEG);
    out = struct('loadedFile', fullfile(R.filepath, R.filename));
    logPrint(R.logFile, sprintf('[load_set] --- Dataset loaded %s/%s ---', R.filepath, R.filename)); % Added [load_set] prefix


end