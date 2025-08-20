function test_pipeline_with_bad_channels()
%TEST_PIPELINE_WITH_BAD_CHANNELS An end-to-end test for the pipeline with bad channel detection.

    try
        fprintf('=== Starting Pipeline with Bad Channel Detection Test ===\n\n');

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
        eog_labels = {'Fp1', 'Fp2'};

        all_labels = {EEG.chanlocs.labels};
        flat_chan_idx = find(strcmp(all_labels, flat_chan_label));
        noisy_chan_idx = find(strcmp(all_labels, noisy_chan_label));
        eog_indices = find(ismember(all_labels, eog_labels));

        fprintf('Injecting flat line into channel: %s (index %d)\n', flat_chan_label, flat_chan_idx);
        EEG.data(flat_chan_idx, 1000:2000) = 0;

        fprintf('Injecting high-amplitude noise into channel: %s (index %d)\n', noisy_chan_label, noisy_chan_idx);
        noise = (rand(1, EEG.pnts) - 0.5) * 500;
        EEG.data(noisy_chan_idx, :) = EEG.data(noisy_chan_idx, :) + noise;
        fprintf('Artifacts injected.\n\n');

        % --- Define Parameters ---
        fprintf('--- Defining Parameters ---\n');
        params = struct();
        params.outputpath = fullfile(pwd, 'tests', 'temp_output');
        bad_chan_params = {...
            'IdxDetect', setdiff(1:EEG.nbchan, eog_indices), ...
            'Kurtosis', true, 'Kurt_Threshold', 4, ...
            'CleanRaw_Flatline', true, 'Flatline_Sec', 5, ...
            'Action', 'remove' ...
        };

        if ~exist(params.outputpath, 'dir'), mkdir(params.outputpath); end
        fprintf('Parameters defined.\n\n');

        % --- Build and Run Pipeline ---
        fprintf('--- Building and Running Pipeline ---\n');
        p = EEGdojo.Params(params);
        pipe = EEGdojo.Pipeline(EEG, p, fullfile(p.outputpath, 'pipeline_bad_chan_test.log'));

        pipe = pipe.addStep(@EEGdojo.preprocess.remove_bad_channels, bad_chan_params{:});

        pipe = pipe.run();
        fprintf('Pipeline finished successfully.\n\n');

        % --- Verification ---
        fprintf('--- Verifying Results ---\n');
        final_EEG = pipe.EEG;
        original_chans = {EEG.urchanlocs.labels};
        final_chans = {final_EEG.chanlocs.labels};

        assert(~ismember(flat_chan_label, final_chans), 'Test failed: Flat channel was not removed by the pipeline.');
        fprintf('OK: Flat channel removed.\n');
        
        assert(~ismember(noisy_chan_label, final_chans), 'Test failed: Noisy channel was not removed by the pipeline.');
        fprintf('OK: Noisy channel removed.\n');

        assert(all(ismember(eog_labels, final_chans)), 'Test failed: EOG channel was incorrectly removed.');
        fprintf('OK: EOG channels were correctly kept.\n');

        % --- Cleanup ---
        rmdir(params.outputpath, 's');

        fprintf('\n=== Pipeline with Bad Channel Detection Test Passed Successfully! ===\n');

    catch ME
        fprintf('\n*** An error occurred during the test: ***\n');
        fprintf('Error identifier: %s\n', ME.identifier);
        fprintf('Error message: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Name: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).name, ME.stack(i).line);
        end
        fprintf('\n=== Testing failed ===\n');
    end

end