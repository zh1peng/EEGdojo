function logplot_badChannels(EEG, badChannels, logPath, detectionType)
    % LOGPLOT_CHANNEL - Save and log properties of detected bad channels.
    %
    % Inputs:
    %   EEG            - EEGLAB EEG structure.
    %   badChannels    - Array of bad channel indices.
    %   logPath        - Path to save channel property figures.
    %   detectionType  - String indicating the detection method (e.g., 'Spec', 'Kurt').

    %
    % Outputs:
    %   None. Saves figures to the specified logPath.

    % Ensure the logPath exists
    if ~exist(logPath, 'dir')
        mkdir(logPath);
    end

    % Loop through detected bad channels
    if ~isempty(badChannels)
        for i = 1:length(badChannels)
            chanIdx = badChannels(i);
                % Plot channel properties using pop_prop
                pop_prop(EEG, 1, chanIdx, NaN, {'freqrange', [2 40]});
                % Save the figure with detection method and channel index
                saveas(gcf, fullfile(logPath, sprintf('BadChannel_%s_%d_Properties.png', detectionType, chanIdx)));
                % Close the figure to avoid clutter
                close(gcf);
        end
    end
end