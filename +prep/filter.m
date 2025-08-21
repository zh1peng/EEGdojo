function [EEG, out] = filter(EEG, varargin)
% FILTER Applies high-pass and low-pass filters to EEG data.
%
% This function applies FIR (Finite Impulse Response) filters to the EEG data
% using EEGLAB's `pop_eegfiltnew` function. It supports both high-pass and
% low-pass filtering.
%
% Inputs:
%   EEG         - EEGLAB EEG structure.
%   varargin    - Optional parameters:
%     'HPfreq'  - (numeric) High-pass cutoff frequency in Hz. If -1 or 0,
%                 no high-pass filter is applied. Default is -1.
%     'LPfreq'  - (numeric) Low-pass cutoff frequency in Hz. If -1 or 0,
%                 no low-pass filter is applied. Default is -1.
%     'LogFile' - (char) Path to the log file for recording processing
%                 information. Default is ''.
%
% Outputs:
%   EEG         - Modified EEGLAB EEG structure with filtered data.
%   out         - Structure containing output information:
%     .HPfreq   - (numeric) The high-pass frequency used.
%     .LPfreq   - (numeric) The low-pass frequency used.
%
% Examples:
%   % 1. Apply a band-pass filter from 0.5 Hz to 30 Hz:
%   EEG = filter(EEG, 'HPfreq', 0.5, 'LPfreq', 30);
%
%   % 2. Apply only a high-pass filter at 1 Hz:
%   EEG = filter(EEG, 'HPfreq', 1);
%
%   % 3. Usage within a pipeline:
%   %    (Assuming 'p' is a parameter structure containing 'p.logFile')
%   pipe = pipe.addStep(@prep.filter, ...
%       'HPfreq', 0.1, ...
%       'LPfreq', 40, ...
%       'LogFile', p.logFile);
%
% See also: pop_eegfiltnew

    % ----------------- Parse inputs -----------------
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('HPfreq', -1, @isnumeric);
    p.addParameter('LPfreq', -1, @isnumeric);
    p.addParameter('LogFile', '', @ischar);

    p.parse(EEG, varargin{:});
    R = p.Results;

    out = struct(); % Initialize output structure
    out.HPfreq = R.HPfreq;
    out.LPfreq = R.LPfreq;

    logPrint(R.LogFile, '[filter] Starting filtering process.');
    
    % Construct log message for filter parameters
    filter_msg = '';
    if R.HPfreq > 0, filter_msg = [filter_msg, sprintf('High-pass: %.2f Hz. ', R.HPfreq)]; end
    if R.LPfreq > 0, filter_msg = [filter_msg, sprintf('Low-pass: %.2f Hz.', R.LPfreq)]; end
    if isempty(filter_msg), filter_msg = 'No filter applied (HPfreq and LPfreq are <= 0).'; end
    logPrint(R.LogFile, sprintf('[filter] %s', filter_msg));

    % Input validation for frequencies
    if R.HPfreq > 0 && R.LPfreq > 0 && R.HPfreq >= R.LPfreq
        error('[filter] High-pass frequency (%.2f Hz) must be lower than low-pass frequency (%.2f Hz).', R.HPfreq, R.LPfreq);
    end
    if R.HPfreq < 0 && R.LPfreq < 0
        logPrint(R.LogFile, '[filter] No valid filter frequencies provided. Skipping filtering.');
        return; % Skip filtering if no valid frequencies
    end

    try
        % Apply high-pass filter if specified
        if R.HPfreq > 0
            EEG = pop_eegfiltnew(EEG, 'locutoff', R.HPfreq, 'plotfreqz', 0);
            EEG = eeg_checkset(EEG);
        end
        
        % Apply low-pass filter if specified
        if R.LPfreq > 0
            EEG = pop_eegfiltnew(EEG, 'hicutoff', R.LPfreq,  'plotfreqz', 0);
            EEG = eeg_checkset(EEG);
        end
        
        logPrint(R.LogFile, '[filter] Filtering complete.');
    catch ME
        error('[filter] Filtering failed: %s', ME.message);
    end
end
