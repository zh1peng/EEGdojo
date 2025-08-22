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
            if isfield(p, 'LogFile'), obj.LogFile = p.LogFile; else, obj.LogFile = 'pipeline.log'; end
            if isfield(p, 'error_LogFile'), obj.error_LogFile = p.error_LogFile; else, obj.error_LogFile = 'pipeline_error.log'; end

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

    function obj = addStep(obj, varargin)
        % ADDSTEP Add a processing step to the pipeline.
        % Usage patterns:
        %   addStep('Name', @func, structArg, 'k1',v1, structArg2, 'k2',v2, ...)
        %   addStep(@func, structArg, 'k1',v1, ...)
        %
        % - Accepts any mix of structs and name-value pairs.
        % - If the first argument is a string/char and the second is a function handle,
        %   that string is used as the step name. Otherwise, name = func2str(func).
        %
        % Example:
        % 1) Named + struct
        % eeg_pipe.addStep('Load', @IO.load_set, localParams.IO);

        % % 2) Unnamed + NV
        % eeg_pipe.addStep(@prep.downsample, 'freq', localParams.DownsamplingRate);

        % % 3) Mix struct + NV (NV overrides struct if same key appears later)
        % eeg_pipe.addStep('Bandpass', @prep.filter, localParams.Filter, 'HighCutoff', 45);

        % % 4) Multiple structs
        % eeg_pipe.addStep('Bad channels', @prep.remove_bad_channels, localParams.BadChan, extraBadChanOpts);

        % % 5) Your original pattern still works
        % eeg_pipe.addStep('Save', @IO.save_set, 'filepath', localParams.IO.outputpath, 'filename', localParams.IO.outfilename);


        assert(~isempty(varargin), 'addStep requires at least a function handle.');

        % Parse leading "name" and func handle
        if (numel(varargin) >= 2) && (ischar(varargin{1}) || isstring(varargin{1})) && isa(varargin{2}, 'function_handle')
            stepName = char(varargin{1});
            func     = varargin{2};
            rest     = varargin(3:end);
        elseif isa(varargin{1}, 'function_handle')
            func     = varargin{1};
            stepName = func2str(func);
            rest     = varargin(2:end);
        else
            error('addStep:InvalidSignature', ...
                'Expected (''Name'', @func, ...) or (@func, ...).');
        end

        % Flatten any combination of structs + NV pairs in "rest"
        args = obj.local_flatten_structs_and_nv(rest);

        % Validate NV pairs
        if mod(numel(args),2) ~= 0
            error('addStep:NameValueMismatch', ...
                  'Arguments after the function must form name–value pairs (got odd count).');
        end
        % Optional: enforce that all names are char/string
        for i = 1:2:numel(args)
            if ~(ischar(args{i}) || isstring(args{i}))
                error('addStep:InvalidKey', 'Name–value key at position %d is not text.', i);
            end
            args{i} = char(args{i}); % normalize to char
        end

        % Store
        s = struct('name', stepName, 'func', func, 'args', {args});
        obj.steps(end+1) = s;
    end


        function obj = addIf(obj, cond, varargin)
            if ~cond, return; end
            obj.addStep(varargin{:});
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
                    logPrint(obj.LogFile, sprintf('[Step %s]', stepName));
                    tStep = tic;

                    % Convention: step functions accept (EEG, p, varargin{:}) and return EEG
                    obj.EEG = step.func(obj.EEG, step.args{:});

                    dt = toc(tStep);
                    logPrint(obj.LogFile, sprintf('[Step %s] completed in %.2fs', stepName, dt));
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



methods (Access = private)
    function flat = local_flatten_structs_and_nv(~, parts)
        % Expand a cell array "parts" that can contain:
        % - structs (expanded to fieldname,value,...),
        % - name–value pairs,
        % - empty [] (ignored).
        % Returns a 1x(2N) cell row of NV pairs.
        flat = {};
        for k = 1:numel(parts)
            a = parts{k};
            if isempty(a)
                continue
            elseif isstruct(a)
                % Allow scalar struct or struct array (concatenate fields in order)
                if numel(a) == 1
                    f = fieldnames(a);
                    v = struct2cell(a);
                    flat = [flat, reshape([f.'; v.'], 1, [])]; %#ok<AGROW>
                else
                    % For struct arrays, expand each element in sequence
                    for j = 1:numel(a)
                        f = fieldnames(a(j));
                        v = struct2cell(a(j));
                        flat = [flat, reshape([f.'; v.'], 1, [])]; %#ok<AGROW>
                    end
                end
            else
                % Pass through (NV pairs or other literals)
                flat = [flat, {a}]; %#ok<AGROW>
            end
        end
    end
end
end

