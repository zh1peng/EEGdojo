classdef Pipeline
    %PIPELINE A class to define and run a sequence of EEG processing steps.

    properties
        EEG
        p
        logFile
        steps
    end

    methods
        function obj = Pipeline(EEG, p, logFile)
            %PIPELINE Construct an instance of this class
            obj.EEG = EEG;
            obj.p = p;
            obj.logFile = logFile;
            obj.steps = {};
            
            if ~exist(fileparts(logFile), 'dir')
                mkdir(fileparts(logFile));
            end
            
            logPrint(obj.logFile, '=== Starting EEG Preprocessing Pipeline ===');
        end

        function obj = addStep(obj, func, varargin)
            %ADDSTEP Add a processing step to the pipeline.
            obj.steps{end+1} = {func, varargin};
        end

        function obj = run(obj)
            %RUN Execute all the steps in the pipeline.
            t1 = tic;
            for i = 1:length(obj.steps)
                step = obj.steps{i};
                func = step{1};
                args = step{2};
                
                try
                    logPrint(obj.logFile, sprintf('--- Running step %d: %s ---', i, func2str(func)));
                    obj.EEG = func(obj.EEG, args{:});
                    logPrint(obj.logFile, sprintf('--- Step %d completed ---', i));
                catch ME
                    logPrint(obj.logFile, sprintf('*** ERROR in step %d (%s): %s ***', i, func2str(func), ME.message));
                    rethrow(ME);
                end
            end
            t2 = toc(t1);
            logPrint(obj.logFile, sprintf('Data preprocessed successfully. It took %.2f seconds.', t2));
            logPrint(obj.logFile, '=== EEG Preprocessing Complete ===');
        end
    end
end

function logPrint(logFile, message)
    fid = fopen(logFile, 'a');
    fprintf(fid, '[%s] %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), message);
    fclose(fid);
    disp(message);
end