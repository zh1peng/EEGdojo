function EEG = remove_powerline(EEG, varargin)
% remove_powerline  Remove mains/line noise via CleanLine or FIR notch.
% Usage:
%   EEG = remove_powerline(EEG, 'Method','cleanline', 'Freq',50, 'BW',2, 'NHarm',3);
%   EEG = remove_powerline(EEG, 'Method','notch',     'Freq',60, 'BW',3, 'NHarm',2);
%
% Params:
%   Method  - 'cleanline' | 'notch' (default 'cleanline')
%   Freq    - fundamental line frequency in Hz (default 50)
%   BW      - half-bandwidth in Hz for notch (±BW) (default 2)
%   NHarm   - number of harmonics to target (default 3). Only those < Nyquist are used.
%   Verbose - true/false (default false)

    % ---- Parse inputs ----
    p = inputParser;
    p.addRequired('EEG', @isstruct);
    p.addParameter('Method','cleanline', @(s) any(strcmpi(s,{'cleanline','notch'})));
    p.addParameter('Freq', 50, @(x) isnumeric(x) && isscalar(x) && x>0);
    p.addParameter('BW',   2,  @(x) isnumeric(x) && isscalar(x) && x>0);
    p.addParameter('NHarm',3,  @(x) isnumeric(x) && isscalar(x) && x>=1);
    p.addParameter('logFile', '', @ischar);
    p.parse(EEG, varargin{:});
    R = p.Results;


    logPrint(R.logFile, sprintf('--- Removing powerline noise using %s method ---', R.Method));

    fs = EEG.srate;
    nyq = fs/2;

    % Build harmonic list under Nyquist
    harm = (1:R.NHarm) * R.Freq;
    harm = harm(harm < nyq);
    if isempty(harm)
        error('[remove_line_noise] No harmonics < Nyquist. \n'); 
    end

    switch lower(R.Method)
        case 'cleanline'
            % --- Adaptive removal using CleanLine ---
            % pop_cleanline accepts vector of line freqs

            logPrint(R.logFile,sprintf('[remove_line_noise] CleanLine at Hz: %s\n', num2str(harm)));
            try
                EEG = pop_cleanline(EEG, 'linefreqs', harm, 'newversion', 1);
                EEG = eeg_checkset(EEG);
                logPrint(R.logFile, '--- CleanLine complete ---');
            catch ME
                error('[remove_line_noise] Error using CleanLine: %s\n', ME.message);
            end

        case 'notch'
            % --- Fixed FIR band-stop around each harmonic: [f-BW, f+BW] ---

            logPrint(R.logFile,'[remove_line_noise] FIR notch (±%.2f Hz) at Hz: %s\n', R.BW, num2str(harm));
            try
                for f0 = harm
                    lo = max(f0 - R.BW, 0);   % lower edge
                    hi = min(f0 + R.BW, nyq); % upper edge
                    if lo <= 0 || hi <= 0 || lo >= hi
                        continue; % skip malformed bands
                    end
                    % pop_eegfiltnew: band-stop when 'revfilt'=1
                    % EEG = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder, revfilt, usefft, plotfreqz)
                    EEG = pop_eegfiltnew(EEG, lo, hi, [], 1, [], 0);
                    EEG = eeg_checkset(EEG);
                end
                logPrint(R.logFile, '--- FIR notch complete ---');
            catch ME
                error('[remove_line_noise] Error using FIR notch: %s\n', ME.message);
            end      
    end
    logPrint(R.logFile, '--- Powerline noise removal complete ---');
end
