classdef Pipeline < handle
    %PIPELINE Define and run a sequence of EEG processing steps (fluent API)

    properties
        EEG
        p
        LogFile
        error_LogFile
        steps % struct array with fields: name, func, args
        continueOnError logical = false
    end

    methods
        function obj = Pipeline(EEG, p)
            if nargin < 1 || isempty(EEG), EEG = []; end
            if nargin < 2 || isempty(p),   p   = struct(); end

            obj.EEG = EEG;
            obj.p   = p;

            % Robust defaults if not provided
            if isfield(p, 'IO') && isfield(p.IO, 'LogFile'),       obj.LogFile = p.IO.LogFile;       else, obj.LogFile = 'pipeline.log'; end
            if isfield(p, 'IO') && isfield(p.IO, 'error_LogFile'), obj.error_LogFile = p.IO.error_LogFile; else, obj.error_LogFile = 'pipeline_error.log'; end

            obj.steps = struct('name', {}, 'func', {}, 'args', {});
            logPrint(obj.LogFile, '=== Starting EEG Preprocessing Pipeline ===');
        end

        function obj = setParam(obj, p)
            obj.p = p;
        end

        function obj = setEEG(obj, EEG)
            obj.EEG = EEG;
        end

        function obj = setIO(obj, LogFile, errLogFile)
            if nargin >= 2 && ~isempty(LogFile),   obj.LogFile = LogFile;     end
            if nargin >= 3 && ~isempty(errLogFile), obj.error_LogFile = errLogFile; end
        end

        function obj = setContinueOnError(obj, tf)
            obj.continueOnError = logical(tf);
        end

        function obj = addStep(obj, func, varargin)
            % addStep(@prep.filter, 'hp',0.1,'lp',40)
            % Optional name: addStep("Filter", @prep.filter, 'hp',0.1)
            if nargin >= 3 && (ischar(func) || isstring(func))
                % User passed a name as first arg; shift
                name = char(func);
                func = varargin{1};
                args = varargin(2:end);
            else
                % Derive a readable name
                name = func2str(func);
                args = varargin;
            end
            s = struct('name', name, 'func', func, 'args', {args});
            obj.steps(end+1) = s;
        end

        function obj = addIf(obj, cond, func, varargin)
            % Conditionally add a step
            if cond
                obj.addStep(func, varargin{:});
            end
        end

        function obj = run(obj, varargin)
            % RUN Execute steps. Options:
            %   'continueOnError', true/false (default: obj.continueOnError)
            ip = inputParser;
            ip.addParameter('continueOnError', obj.continueOnError, @(x)islogical(x)&&isscalar(x));
            ip.parse(varargin{:});
            contErr = ip.Results.continueOnError;

            t1 = tic;
            for i = 1:numel(obj.steps)
                step = obj.steps(i);
                stepName = sprintf('#%d %s', i, step.name);
                try
                    logPrint(obj.LogFile, sprintf('--- Running %s ---', stepName));
                    tStep = tic;

                    % Convention: step functions accept (EEG, p, varargin{:}) and return EEG
                    obj.EEG = step.func(obj.EEG, step.args{:});

                    dt = toc(tStep);
                    logPrint(obj.LogFile, sprintf('--- %s completed in %.2fs ---', stepName, dt));
                catch ME
                    logPrint(obj.error_LogFile, sprintf('*** ERROR in %s: %s ***', stepName, ME.message));
                    if contErr
                        % Keep going; mark in log
                        logPrint(obj.LogFile, sprintf('!!! Continuing after error in %s', stepName));
                    else
                        rethrow(ME);
                    end
                end
            end
            t2 = toc(t1);
            logPrint(obj.LogFile, sprintf('Data preprocessed successfully. Total time: %.2f seconds.', t2));
            logPrint(obj.LogFile, '=== EEG Preprocessing Complete ===');
        end
    end
end
