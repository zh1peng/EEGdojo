function IO = setIOPaths(IO)
    % setIOPaths - Set input/output paths for EEG preprocessing
    %
    % Syntax: IO = setIOPaths(IO)
    %
    % Inputs:
    %    IO - Structure containing input/output paths and filenames
    %
    % Outputs:
    %    IO - Updated structure with additional paths and filenames

    % Extract the base name from the input filename
    if isempty(IO.filename)
        print('No input filename specified.');
        return;
    end

    [~, baseName, ~] = fileparts(IO.filename);

    % Define the log path and output filenames
    IO.logPath = fullfile(IO.outputpath, [baseName, '_log']);
    IO.outfilename = [baseName '_preprocessed.set'];
    IO.ECGfilename = [baseName '_ECG.set'];
    IO.error_LogFile = fullfile(IO.logPath, [baseName, '_error.log']);
    IO.LogFile = fullfile(IO.logPath, [baseName '_preprocessed.log']);

    % Delete existing log files if they exist
    if exist(IO.error_LogFile, 'file')
        delete(IO.error_LogFile);
    end
    if exist(IO.LogFile, 'file')
        delete(IO.LogFile);
    end
    
    IO.baseName = baseName;

    % Check if the output path exists, if not, create it
    if ~exist(IO.outputpath, 'dir')
        status = mkdir(IO.outputpath);
        if ~status
            error('Failed to create output directory: %s', IO.outputpath);
        end
    end

    % Check if the log path exists, if not, create it
    if ~exist(IO.logPath, 'dir')
        status = mkdir(IO.logPath);
        if ~status
            error('Failed to create log directory: %s', IO.logPath);
        end
    end
end