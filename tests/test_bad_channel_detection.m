function test_bad_channel_detection()
%TEST_BAD_CHANNEL_DETECTION A specific test for the bad channel detection utility.

    try
        fprintf('=== Starting Bad Channel Detection Test ===\n\n');

        % --- Setup ---
        fprintf('--- Setup ---\n');
        addpath(fileparts(pwd));
        test_data_path = 'sub-101_task-moviep_run-01_eeg.set';
        fprintf('Loading test data: %s\n', test_data_path);
        EEG = pop_loadset('filename', test_data_path, 'filepath', pwd);
        fprintf('Test data loaded successfully.\n\n');

        % --- Artificially Create Bad Channels ---
        fprintf('--- Creating Artificial Artifacts ---\n');
        
        flat_chan_label = 'C3';
        noisy_chan_label = 'Pz';
        known_bad_label = 'O1';
        eog_labels = {'Fp1', 'Fp2'};

        all_labels = {EEG.chanlocs.labels};
        flat_chan_idx = find(strcmp(all_labels, flat_chan_label));
        noisy_chan_idx = find(strcmp(all_labels, noisy_chan_label));
        known_bad_idx = find(strcmp(all_labels, known_bad_label));
        eog_indices = find(ismember(all_labels, eog_labels));

        fprintf('Injecting flat line into channel: %s (index %d)\n', flat_chan_label, flat_chan_idx);
        EEG.data(flat_chan_idx, 1000:2000) = 0;

        fprintf('Injecting high-amplitude noise into channel: %s (index %d)\n', noisy_chan_label, noisy_chan_idx);
        noise = (rand(1, EEG.pnts) - 0.5) * 500;
        EEG.data(noisy_chan_idx, :) = EEG.data(noisy_chan_idx, :) + noise;
        fprintf('Artifacts injected.\n\n');

        % --- Run Bad Channel Detection ---
        fprintf('--- Running Bad Channel Detection ---\n');
        
        [~, out] = EEGdojo.preprocess.remove_bad_channels(EEG, ...
            'IdxDetect', setdiff(1:EEG.nbchan, eog_indices), ...
            'Kurtosis', true, 'Kurt_Threshold', 4, ...
            'CleanRaw_Flatline', true, 'Flatline_Sec', 5, ...
            'KnownBadIdx', known_bad_idx, ...
            'Action', 'stash');

        fprintf('Detection complete. Found %d bad channels.\n\n', length(out.Bad.all));

        % --- Verification ---
        fprintf('--- Verifying Results ---\n');
        
        assert(ismember(flat_chan_idx, out.Bad.all), 'Test failed: Flat channel was not detected.');
        fprintf('OK: Flat channel detected.\n');
        
        assert(ismember(noisy_chan_idx, out.Bad.all), 'Test failed: Noisy channel was not detected.');
        fprintf('OK: Noisy channel detected.\n');

        assert(ismember(known_bad_idx, out.Bad.all), 'Test failed: Known bad channel was not included.');
        fprintf('OK: Known bad channel included.\n');

        assert(~any(ismember(eog_indices, out.Bad.all)), 'Test failed: EOG channel was incorrectly flagged as bad.');
        fprintf('OK: EOG channels were correctly ignored.\n');

        fprintf('\n=== Bad Channel Detection Test Passed Successfully! ===\n');

    catch ME
        fprintf('\n*** An error occurred during the bad channel detection test: ***\n');
        fprintf('Error identifier: %s\n', ME.identifier);
        fprintf('Error message: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Name: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).name, ME.stack(i).line);
        end
        fprintf('\n=== Testing failed ===\n');
    end

end