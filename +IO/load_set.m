function [EEG, out] = load_set(~, varargin)
% LOAD_SET  Load an EEGLAB .set dataset
%
% Usage:
%   [EEG, out] = eegdojo.io.load_set([], 'filename','x.set','filepath','./data')
%
    p = inputParser;
    p.addParameter('filename','',@ischar);
    p.addRequired(p,'filepath','',@ischar);
    p.addRequired(p,'logFile', '', @ischar);
    p.parse(varargin{:});

    R = p.Results;
    logPrint(R.logFile, sprintf('--- Loading dataset: %s/%s ---', R.filepath, R.filename));
    try
        EEG = pop_loadset('filename', R.filename, 'filepath', R.filepath);
        EEG = eeg_checkset(EEG);
        out = struct('loadedFile', fullfile(R.filepath, R.filename));
        logPrint(R.logFile, sprintf('--- Dataset loaded %s/%s ---', R.filepath, R.filename));
    catch ME
       error('[load_set] Error loading dataset: %s', ME.message)
    end

    
end
