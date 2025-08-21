function run_all_tests_simple()
%RUN_ALL_TESTS_SIMPLE A simple script to test the EEGdojo toolbox.

    try
        fprintf('=== Starting EEGdojo Test Script ===\n\n');

        % --- Setup ---
        fprintf('--- Setup ---\n');
        % Add project root to path
        % Assuming EEGdojo root is the parent directory of 'tests'
        eegdojo_root = fileparts(pwd);
        addpath(genpath(eegdojo_root));

        % Add EEGLAB and FASTER paths (adjust these paths as necessary for your system)
        % For testing, you might need to manually set these or ensure they are in MATLAB's path
        addpath(genpath('Z:\matlab_toolbox\eeglab2023.1'));
        addpath(genpath('Z:\matlab_toolbox\FASTER'));

        test_data_path = fullfile(eegdojo_root, 'tests', 'data', 'test_data.set');
        fprintf('Loading test data: %s\n', test_data_path);
        % Ensure the test data exists. If not, you might need to provide a dummy EEG struct.
        if ~exist(test_data_path, 'file')
            warning('Test data not found at %s. Creating a dummy EEG structure for testing.', test_data_path);
            EEG = pop_EEG_dummy(); % A dummy function to create a basic EEG struct
        else
            EEG = pop_loadset('filename', 'test_data.set', 'filepath', fullfile(eegdojo_root, 'tests', 'data'));
        end
        original_EEG = EEG; % Keep a copy for tests that modify EEG
        fprintf('Test data loaded/created successfully.\n\n');
        
        original_EEG=prep.downsample(original_EEG, 'freq', 150);

        % --- Test +IO functions ---
        fprintf('--- Testing +IO functions ---\n');

        % Test filesearch_regexp
        disp('Testing filesearch_regexp...');
        files = filesearch_regexp(fullfile(eegdojo_root, 'tests', 'data'), '.*\.set');
        assert(~isempty(files), 'filesearch_regexp failed: no files found.');
        assert(iscell(files), 'filesearch_regexp failed: did not return a cell array.');
        disp('OK');

        % Test save_set and load_set
        disp('Testing IO.save_set and IO.load_set...');
        save_filename = 'temp_test_save.set';
        temp_dir = fullfile(eegdojo_root, 'tests', 'temp');
        if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end
        IO.save_set(EEG, 'filename', save_filename, 'filepath', temp_dir);
        assert(exist(fullfile(temp_dir, save_filename), 'file') == 2, 'save_set failed: file not created.');
        [loaded_EEG, ~] = IO.load_set([], 'filename', save_filename, 'filepath', temp_dir);
        assert(isequal(EEG.data, loaded_EEG.data), 'load_set failed: loaded data does not match original data.');
        delete(fullfile(temp_dir, save_filename));
        rmdir(temp_dir); % Clean up temp directory
        disp('OK');
        fprintf('\n');

        % --- Test +prep functions ---
        fprintf('--- Testing +prep functions ---\n');

        % Test downsample
        disp('Testing prep.downsample...');
        EEG_ds = prep.downsample(original_EEG, 'freq', 125); % Use original_EEG
        assert(EEG_ds.srate == 125, 'downsample failed: incorrect srate.');
        disp('OK');

      
        % Test filter
        disp('Testing prep.filter...');
        EEG_filt = prep.filter(original_EEG, 'LowCutoff', 1, 'HighCutoff', 30); % Use original_EEG
        assert(~isempty(EEG_filt.data), 'filter failed: data is empty.');
        disp('OK');

        % Test remove_channels (updated)
        disp('Testing prep.remove_channels...');
        % Test by labels
        chans_to_remove_labels = {'Fp1', 'Fp2'};
        iscellstr(chans_to_remove_labels)
        EEG_rem_labels = prep.remove_channels(original_EEG, 'ChanLabels', chans_to_remove_labels);
        assert(~ismember('Fp1', {EEG_rem_labels.chanlocs.labels}), 'remove_channels by labels failed.');
        % Test by indices
        chans_to_remove_idx = [1 2]; % Assuming Fp1 and Fp2 are first two channels
        EEG_rem_idx = prep.remove_channels(original_EEG, 'ChanIdx', chans_to_remove_idx);
        assert(EEG_rem_idx.nbchan == original_EEG.nbchan - length(chans_to_remove_idx), 'remove_channels by index failed.');
        disp('OK');

        % Test reref
        disp('Testing prep.reref...');
        EEG_reref = prep.reref(original_EEG); % Use original_EEG
        assert(size(EEG_reref.data,1) == size(original_EEG.data,1), 'reref failed.');
        disp('OK');

        % Test remove_powerline
        disp('Testing prep.remove_powerline...');
        EEG_pl = prep.remove_powerline(original_EEG, 'Freq', 50); % Use original_EEG
        assert(~isempty(EEG_pl.data), 'remove_powerline failed.');
        disp('OK');

        % Test interpolate
        disp('Testing prep.interpolate...');
        % Create a dummy EEG with a missing channel for interpolation test
        EEG_interp_test = original_EEG;
        if EEG_interp_test.nbchan > 1
            missing_chan_label = EEG_interp_test.chanlocs(1).labels;
            EEG_interp_test = pop_select(EEG_interp_test, 'nochannel', 1);
            EEG_interp_test.urchanlocs = original_EEG.chanlocs; % Set urchanlocs for interpolate
            EEG_interp = prep.interpolate(EEG_interp_test);
            assert(ismember(missing_chan_label, {EEG_interp.chanlocs.labels}), 'interpolate failed: channel not re-added.');
        else
            warning('Skipping interpolate test: Not enough channels in dummy EEG.');
        end
        disp('OK');

        % Test remove_bad_channels
        disp('Testing prep.remove_bad_channels...');
        % This test is tricky as it depends on data characteristics.
        % We'll enable a simple detector and check if channel count changes or mask is created.
        EEG_badchan = prep.remove_bad_channels(original_EEG, 'Kurtosis', true, 'Kurt_Threshold', 3, 'Action', 'remove');
        assert(EEG_badchan.nbchan <= original_EEG.nbchan, 'remove_bad_channels failed: channel count increased.');
        EEG_badchan_flag = prep.remove_bad_channels(original_EEG, 'Kurtosis', true, 'Kurt_Threshold', 3, 'Action', 'flag');
        assert(isfield(EEG_badchan_flag.etc, 'clean_channel_mask'), 'remove_bad_channels failed: clean_channel_mask not created.');
        disp('OK');

        % Test remove_bad_ICs
        disp('Testing prep.remove_bad_ICs...');
        % First, ensure ICA is computed or will be computed by the function
        EEG_ica_test = original_EEG;
        if ~isfield(EEG_ica_test, 'icaweights') || isempty(EEG_ica_test.icaweights)
            % If ICA not present, run a simple ICA
            EEG_ica_test = pop_runica(EEG_ica_test, 'extended', 1, 'pca', min(EEG_ica_test.nbchan-1, 30));
        end
        % Add a dummy ECG channel for testing DetectECG
        if ~ismember('ECG', {EEG_ica_test.chanlocs.type})
            % Add a dummy ECG channel if not present
            ecg_chan = EEG_ica_test.chanlocs(1); % Copy first channel info
            ecg_chan.labels = 'ECG';
            ecg_chan.type = 'ECG';
            EEG_ica_test.chanlocs(end+1) = ecg_chan;
            EEG_ica_test.data(end+1,:) = randn(1, EEG_ica_test.pnts); % Dummy ECG data
            EEG_ica_test.nbchan = EEG_ica_test.nbchan + 1;
            EEG_ica_test = eeg_checkset(EEG_ica_test);
        end

        EEG_badic = prep.remove_bad_ICs(EEG_ica_test, 'ICLabelOn', true, 'FASTEROn', true, 'DetectECG', true); % Explicitly set DetectECG to true for test
        assert(isfield(EEG_badic, 'icaweights'), 'remove_bad_ICs failed: ICA weights missing.');
        assert(isfield(EEG_badic.etc, 'EEGdojo'), 'remove_bad_ICs failed: EEGdojo field missing.');
        disp('OK');

        % Test segment_rest
        disp('Testing prep.segment_rest...');
        EEG_rest_seg = prep.segment_rest(original_EEG, 'EpochLength', 1, 'EpochOverlap', 0.5);
        assert(EEG_rest_seg.trials > 0, 'segment_rest failed: no epochs created.');
        assert(isfield(EEG_rest_seg, 'epoch'), 'segment_rest failed: epoch field missing.');
        disp('OK');

        % Test segment_task
        disp('Testing prep.segment_task...');
        % Add dummy events if none exist or if specific markers are needed
        if isempty(original_EEG.event)
            EEG_task_test = original_EEG;
            % Add some dummy events
            EEG_task_test.event(1).type = 'stim_on';
            EEG_task_test.event(1).latency = 1000;
            EEG_task_test.event(2).type = 'resp';
            EEG_task_test.event(2).latency = 2000;
            EEG_task_test.event(3).type = 'stim_on';
            EEG_task_test.event(3).latency = 3000;
            EEG_task_test = eeg_checkset(EEG_task_test);
        else
            EEG_task_test = original_EEG;
        end
        
        % Ensure there's at least one event of type 'stim_on' or 'resp'
        if ~ismember('stim_on', {EEG_task_test.event.type}) && ~ismember('resp', {EEG_task_test.event.type})
            % If no relevant events, add one
            EEG_task_test.event(end+1).type = 'stim_on';
            EEG_task_test.event(end).latency = 1000;
            EEG_task_test = eeg_checkset(EEG_task_test);
        end

        EEG_task_seg = prep.segment_task(EEG_task_test, 'Markers', {'stim_on', 'resp'}, 'TimeWindow', [-0.1 0.5]);
        assert(EEG_task_seg.trials > 0, 'segment_task failed: no epochs created.');
        assert(isfield(EEG_task_seg, 'epoch'), 'segment_task failed: epoch field missing.');
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

% Helper function to create a dummy EEG structure if test data is not found
function EEG = pop_EEG_dummy()
    EEG = eeg_emptyset();
    EEG.srate = 250;
    EEG.nbchan = 32;
    EEG.pnts = 250 * 10; % 10 seconds of data
    EEG.trials = 1;
    EEG.data = randn(EEG.nbchan, EEG.pnts);
    EEG.chanlocs = struct('labels', cell(1, EEG.nbchan), 'type', 'EEG');
    for i = 1:EEG.nbchan
        EEG.chanlocs(i).labels = sprintf('Chan%d', i);
    end
    EEG = eeg_checkset(EEG);
end