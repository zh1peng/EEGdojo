function [EEG, out] = save_set(EEG, varargin)
% SAVE_SET  Save EEGLAB dataset (.set)
%
% Usage:
%   [EEG, out] = eegdojo.io.save_set(EEG, 'filename','x.set','filepath','./out')
%
    p = inputParser;
    p.addRequired(p,'EEG',@isstruct);
    p.addRequired(p,'filename','',@ischar);
    p.addRequired(p,'filepath','',@ischar);
    p.addParameter('logFile', '', @ischar);
    p.parse(varargin{:});
    R = p.Results;
     sprintf('--- Saving dataset: %s/%s ---',R.filepath, R.filename);
    try
        pop_saveset(EEG, 'filename', R.filename, 'filepath', R.filepath);
        out = struct('savedFile', fullfile(R.filepath, R.filename));
    logPrint(R.logFile, sprintf('--- Dataset saved: %s/%s ---', R.filepath, R.filename));
    catch ME
        error('[save_set] Error saving dataset: %s', ME.message)
    end
    
end
