function [EEG, out] = save_set(EEG, varargin)
% SAVE_SET  Save EEGLAB dataset (.set)
%
% Usage:
%   [EEG, out] = eegdojo.io.save_set(EEG, 'filename','x.set','filepath','./out')
%
    p = inputParser;
    p.addParameter('filename','',@ischar);
    p.addParameter('filepath','',@ischar);
    p.addParameter('logFile', '', @ischar);
    p.addParameter('error_logFile', '', @ischar);
    p.parse(varargin{:});
    R = p.Results;
    logPrint(R.logFile, sprintf('--- Saving dataset: %s ---', R.filename));
    try
        pop_saveset(EEG, 'filename', R.filename, 'filepath', R.filepath);
        out = struct('savedFile', fullfile(R.filepath, R.filename));
    catch ME
        logPrint(R.error_logFile, '[save_set] Error saving dataset: %s', ME.message)
        return;
    end
    
end
