function test_pipeline()
%TEST_PIPELINE An end-to-end test for the EEGdojo pipeline.

    try
        fprintf('=== Starting EEGdojo Pipeline Test ===\n\n');

        % --- Setup ---
        fprintf('--- Setup ---\n');
        addpath(fileparts(pwd));
        test_data_path = 'sub-101_task-moviep_run-01_eeg.set';
        fprintf('Loading test data: %s\n', test_data_path);
        EEG = pop_loadset('filename', test_data_path, 'filepath', pwd);
        fprintf('Test data loaded successfully.\n\n');

        % --- Define Parameters ---
        fprintf('--- Defining Parameters ---\n');
        params = struct();
        params.outputpath = fullfile(pwd, 'tests', 'temp_output');
        params.filename = test_data_path;
        params.filepath = pwd;
        params.chan2remove = {'Fp1', 'Fp2'};
        params.DownsamplingRate = 128;
        params.HPfreq = 1;
        params.LPfreq = 40;
        params.PowerLineFrequency = 50;
        % For this test, we will skip the more complex steps
        params.StartMarker = '1';
        params.EndMarker = '2';
        params.eogchan = {};
        params.ICLabel = [];

        % Create output directory
        if ~exist(params.outputpath, 'dir')
            mkdir(params.outputpath);
        end
        fprintf('Parameters defined.\n\n');

        % --- Build and Run Pipeline ---
        fprintf('--- Building and Running Pipeline ---\n');
        p = EEGdojo.Params(params);
        pipe = EEGdojo.Pipeline(EEG, p, fullfile(p.outputpath, 'pipeline_test.log'));

        pipe = pipe.addStep(@EEGdojo.preprocess.remove_channels, 'labels', p.chan2remove);
        pipe = pipe.addStep(@EEGdojo.preprocess.downsample, 'freq', p.DownsamplingRate);
        pipe = pipe.addStep(@EEGdojo.preprocess.filter, 'HPfreq', p.HPfreq, 'LPfreq', p.LPfreq);
        pipe = pipe.addStep(@EEGdojo.preprocess.remove_powerline, 'Freq', p.PowerLineFrequency);
        pipe = pipe.addStep(@EEGdojo.preprocess.reref);

        pipe = pipe.run();
        fprintf('Pipeline finished successfully.\n\n');

        % --- Verification ---
        fprintf('--- Verifying Results ---\n');
        final_EEG = pipe.EEG;
        assert(final_EEG.srate == p.DownsamplingRate, 'Pipeline test failed: incorrect sampling rate.');
        assert(~ismember('Fp1', {final_EEG.chanlocs.labels}), 'Pipeline test failed: channel removal did not work.');
        fprintf('All pipeline checks passed!\n\n');

        % --- Cleanup ---
        rmdir(params.outputpath, 's');

        fprintf('=== Pipeline Test Complete ===\n');

    catch ME
        fprintf('\n*** An error occurred during the pipeline test: ***\n');
        fprintf('Error identifier: %s\n', ME.identifier);
        fprintf('Error message: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Name: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).name, ME.stack(i).line);
        end
        fprintf('\n=== Pipeline testing failed ===\n');
    end

end
