function run_all_tests_simple()
%RUN_ALL_TESTS_SIMPLE A simple script to test the EEGdojo toolbox.

    try
        fprintf('=== Starting EEGdojo Test Script ===\n\n');

        % --- Setup ---
        fprintf('--- Setup ---\n');
        % Add project root to path
        addpath(genpath('Z:\matlab_toolbox\eeglab2023.1'))
        addpath(genpath('Z:\matlab_toolbox\FASTER'))
        addpath(genpath(fileparts(pwd)));
        test_data_path = 'sub-101_task-moviep_run-01_eeg.set';
        fprintf('Loading test data: %s\n', test_data_path);
        EEG = pop_loadset('filename', test_data_path, 'filepath', pwd);
        fprintf('Test data loaded successfully.\n\n');

        % --- Test +IO functions ---
        fprintf('--- Testing +IO functions ---\n');
        
        % Test filesearch_regexp
       
        files = filesearch_regexp(pwd, '.*\.set');
        assert(~isempty(files), 'filesearch_regexp failed: no files found.');
        assert(iscell(files), 'filesearch_regexp failed: did not return a cell array.');
        disp('OK');

        % Test save_set and load_set
        disp('Testing EEGdojo.IO.save_set and EEGdojo.IO.load_set...');
        save_filename = 'temp_test_save.set';
        EEGdojo.IO.save_set(EEG, 'filename', save_filename, 'filepath', pwd);
        assert(exist(fullfile(pwd, save_filename), 'file') == 2, 'save_set failed: file not created.');
        [loaded_EEG, ~] = EEGdojo.IO.load_set([], 'filename', save_filename, 'filepath', pwd);
        assert(isequal(EEG.data, loaded_EEG.data), 'load_set failed: loaded data does not match original data.');
        delete(fullfile(pwd, save_filename));
        disp('OK');
        fprintf('\n');

        % --- Test +preprocess functions ---
        fprintf('--- Testing +preprocess functions ---\n');

        % Test downsample
        disp('Testing EEGdojo.preprocess.downsample...');
        EEG_ds = EEGdojo.preprocess.downsample(EEG, 'freq', 125);
        assert(EEG_ds.srate == 125, 'downsample failed: incorrect srate.');
        disp('OK');

        % Test filter
        disp('Testing EEGdojo.preprocess.filter...');
        EEG_filt = EEGdojo.preprocess.filter(EEG, 'HPfreq', 1, 'LPfreq', 30);
        assert(~isempty(EEG_filt.data), 'filter failed: data is empty.');
        disp('OK');

        % Test remove_channels
        disp('Testing EEGdojo.preprocess.remove_channels...');
        chans_to_remove = {'Fp1', 'Fp2'};
        EEG_rem = EEGdojo.preprocess.remove_channels(EEG, 'labels', chans_to_remove);
        assert(~ismember('Fp1', {EEG_rem.chanlocs.labels}), 'remove_channels failed.');
        disp('OK');

        % Test reref
        disp('Testing EEGdojo.preprocess.reref...');
        EEG_reref = EEGdojo.preprocess.reref(EEG);
        assert(size(EEG_reref.data,1) == size(EEG.data,1), 'reref failed.');
        disp('OK');

        % Test remove_powerline
        disp('Testing EEGdojo.preprocess.remove_powerline...');
        EEG_pl = EEGdojo.preprocess.remove_powerline(EEG, 'Freq', 50);
        assert(~isempty(EEG_pl.data), 'remove_powerline failed.');
        disp('OK');

        % Test interpolate
        disp('Testing EEGdojo.preprocess.interpolate...');
        urchanlocs = EEG.chanlocs;
        EEG_to_interp = pop_select(EEG, 'nochannel', {'Fz'});
        EEG_to_interp.urchanlocs = urchanlocs;
        EEG_interp = EEGdojo.preprocess.interpolate(EEG_to_interp);
        assert(ismember('Fz', {EEG_interp.chanlocs.labels}), 'interpolate failed.');
        disp('OK');

        fprintf('\n=== All checks passed successfully! ===\n');

    catch ME
        fprintf('\n*** An error occurred during testing: ***\n');
        fprintf('Error identifier: %s\n', ME.identifier);
        fprintf('Error message: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Name: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).name, ME.stack(i).line);
        end
        fprintf('\n=== Testing failed ===\n');
    end

end
