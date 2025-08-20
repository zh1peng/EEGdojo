function [EEG, out] = load_set(~, varargin)
% LOAD_SET  Load an EEGLAB .set dataset
%
% Usage:
%   [EEG, out] = eegdojo.io.load_set([], 'filename','x.set','filepath','./data')
%
    p = inputParser;
    p.addParameter('filename','',@ischar);
    p.addParameter('filepath','',@ischar);
    p.addParameter('logFile', '', @ischar);
    p.addParameter('error_logFile', '', @ischar);
    p.parse(varargin{:});

    R = p.Results;
    logPrint(R.logFile, sprintf('--- Loading dataset: %s ---', R.filename));
    try
        EEG = pop_loadset('filename', R.filename, 'filepath', R.filepath);
        EEG = eeg_checkset(EEG);
        out = struct('loadedFile', fullfile(R.filepath, R.filename));
        logPrint(R.logFile, '--- Dataset loaded ---');
    catch ME
        logPrint(R.error_logFile, '[load_set] Error loading dataset: %s', ME.message)
        return;
    end

    
end
