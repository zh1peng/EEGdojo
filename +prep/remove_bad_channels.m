function [EEG, out] = remove_bad_channels(EEG, varargin)
% REMOVE_BAD_CHANNELS  Modular bad-channel detection (kurt/spec/prob + FASTER + CleanRaw).
% Pipeline-friendly: returns [EEG,out], with optional immediate removal/flagging.
%
% Example (pipeline):
%   pipe = pipe.addStep(@eegdojo.preprocess.remove_bad_channels, ...
%       'IdxDetect', idxEEG, ...
%       'Kurtosis',true,'Kurt_Threshold',5, ...
%       'Probability',true,'Prob_Threshold',5, ...
%       'Spectrum',true,'Spec_Threshold',5,'Spec_FreqRange',[1 45], ...
%       'FASTER_MeanCorr',true,'FASTER_Threshold',0.4,'FASTER_RefChan',[], ...
%       'FASTER_Variance',true,'FASTER_VarThreshold',3, ...
%       'FASTER_Hurst',true,'FASTER_HurstThreshold',3, ...
%       'CleanRaw_Flatline',true,'Flatline_Sec',5, ...
%       'CleanRaw_Noise',true,'CleanChan_Corr',0.8,'CleanChan_Line',4, ...
%       'KnownBadIdx', knownBad, ...
%       'PlotFcn',@(EEG,idx,logPath,tag) EEGdojo_logplot_badChannels(EEG,idx,logPath,tag), ...
%       'LogPath', p.outputpath, ...
%       'Action','stash');   % 'stash' | 'remove' | 'flag'
%
% Assumptions:
%   - FASTER_rejchan and cleanraw_rejchan are available and accept the used args.
%   - pop_rejchan exists (EEGLAB).


    % --------- Parse inputs ----------
    p = inputParser;
    p.addRequired('EEG', @isstruct);

    % scope / IO / action
    p.addParameter('IdxDetect', [], @(x) isnumeric(x) || islogical(x));
    p.addParameter('Action', 'remove', @(s) any(strcmpi(s,{'remove','flag'})));
    p.addParameter('PlotFcn', [], @(f) isempty(f) || isa(f,'function_handle'));
    p.addParameter('LogPath', pwd, @(s) ischar(s) || isstring(s));
    p.addParameter('KnownBadIdx', [], @(x) isempty(x) || isnumeric(x));

    % classic EEGLAB detectors
    p.addParameter('Kurtosis', false, @islogical);
    p.addParameter('Kurt_Threshold', 5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    p.addParameter('Probability', false, @islogical);
    p.addParameter('Prob_Threshold', 5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    p.addParameter('Spectrum', false, @islogical);
    p.addParameter('Spec_Threshold', 5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    p.addParameter('Spec_FreqRange', [1 50], @(x)isnumeric(x)&&numel(x)==2&&all(x>=0));
    p.addParameter('NormOn', 'on', @(x) any(strcmpi(x,{'on','off'}))); % for pop_rejchan

    % FASTER (MeanCorr / Variance / Hurst)
    p.addParameter('FASTER_MeanCorr', false, @islogical);
    p.addParameter('FASTER_Threshold', 0.4, @isnumeric);     % r-cutoff for MeanCorr
    p.addParameter('FASTER_RefChan', [], @(x) isempty(x) || (isscalar(x)&&isnumeric(x)));
    p.addParameter('FASTER_Bandpass', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
    p.addParameter('FASTER_Variance', false, @islogical);
    p.addParameter('FASTER_VarThreshold', 3, @isnumeric);    % typical z/FASTER cutoff
    p.addParameter('FASTER_Hurst', false, @islogical);
    p.addParameter('FASTER_HurstThreshold', 3, @isnumeric);

    % CleanRaw: flatline
    p.addParameter('CleanRaw_Flatline', false, @islogical);
    p.addParameter('Flatline_Sec', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
    p.addParameter('CleanDrift_Band', [0.25 0.75], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));

    % CleanRaw: noisy channels (corr + line)
    p.addParameter('CleanRaw_Noise', false, @islogical);
    p.addParameter('CleanChan_Corr', 0.8, @isnumeric);
    p.addParameter('CleanChan_Line', 4, @isnumeric);
    p.addParameter('CleanChan_MaxBad', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
    p.addParameter('CleanChan_NSamp', 50, @isnumeric);

    p.parse(EEG, varargin{:});
    R = p.Results;

    if isempty(R.IdxDetect), R.IdxDetect = 1:EEG.nbchan; end
    if islogical(R.IdxDetect), R.IdxDetect = find(R.IdxDetect); end
    if ~exist(R.LogPath,'dir'), mkdir(R.LogPath); end

    % Preserve original channel locs, keep nodatchans field sane
    if ~isfield(EEG,'urchanlocs') || isempty(EEG.urchanlocs), [EEG.urchanlocs] = deal(EEG.chanlocs); end
    if ~isfield(EEG,'chaninfo') || ~isfield(EEG.chaninfo,'nodatchans'), EEG.chaninfo.nodatchans = []; end

    % --------- Run detectors ----------
    Bad  = struct();
    used = {};

    % 1) Kurtosis
    if R.Kurtosis
        [~, idxRel] = pop_rejchan(EEG, 'elec', R.IdxDetect, ...
            'threshold', R.Kurt_Threshold, 'norm', lower(R.NormOn), 'measure', 'kurt');
        Bad.Kurt = R.IdxDetect(idxRel);
        used{end+1} = 'Kurt'; doPlot(R.PlotFcn, EEG, Bad.Kurt, R.LogPath, 'Kurt');
    else
        Bad.Kurt = []; 
    end

    % 2) Spectrum
    if R.Spectrum
        [~, idxRel] = pop_rejchan(EEG, 'elec', R.IdxDetect, ...
            'threshold', R.Spec_Threshold, 'norm', lower(R.NormOn), ...
            'measure', 'spec', 'freqrange', R.Spec_FreqRange);
        Bad.Spec = R.IdxDetect(idxRel);
        used{end+1} = 'Spec'; doPlot(R.PlotFcn, EEG, Bad.Spec, R.LogPath, 'Spec');
    else
        Bad.Spec = []; 
    end

    % 3) Probability
    if R.Probability
        [~, idxRel] = pop_rejchan(EEG, 'elec', R.IdxDetect, ...
            'threshold', R.Prob_Threshold, 'norm', lower(R.NormOn), 'measure', 'prob');
        Bad.Prob = R.IdxDetect(idxRel);
        used{end+1} = 'Prob'; doPlot(R.PlotFcn, EEG, Bad.Prob, R.LogPath, 'Prob');
    else
        Bad.Prob = []; 
    end

    % 4) FASTER: MeanCorr
    if R.FASTER_MeanCorr
        args = {'elec', R.IdxDetect, 'measure','meanCorr', 'threshold', R.FASTER_Threshold};
        if ~isempty(R.FASTER_RefChan), args = [args, {'refchan', R.FASTER_RefChan}]; end
        if ~isempty(R.FASTER_Bandpass), args = [args, {'bandpass', R.FASTER_Bandpass}]; end
        [~, idx] = FASTER_rejchan(EEG, args{:});
        % exclude refchan if included
        if ~isempty(R.FASTER_RefChan)
            idx = setdiff(idx, R.FASTER_RefChan);
        end
        Bad.MeanCorr = idx(:)';
        used{end+1} = 'FASTER_MeanCorr'; doPlot(R.PlotFcn, EEG, Bad.MeanCorr, R.LogPath, 'MeanCorr');
    else
        Bad.MeanCorr = [];
    end

    % 5) FASTER: Variance
    if R.FASTER_Variance
        args = {'elec', R.IdxDetect, 'measure','variance', 'threshold', R.FASTER_VarThreshold};
        if ~isempty(R.FASTER_RefChan), args = [args, {'refchan', R.FASTER_RefChan}]; end
        if ~isempty(R.FASTER_Bandpass), args = [args, {'bandpass', R.FASTER_Bandpass}]; end
        [~, idx] = FASTER_rejchan(EEG, args{:});
        if ~isempty(R.FASTER_RefChan)
            idx = setdiff(idx, R.FASTER_RefChan);
        end
        Bad.Variance = idx(:)';
        used{end+1} = 'FASTER_Variance'; doPlot(R.PlotFcn, EEG, Bad.Variance, R.LogPath, 'Variance');
    else
        Bad.Variance = [];
    end

    % 6) FASTER: Hurst
    if R.FASTER_Hurst
        args = {'elec', R.IdxDetect, 'measure','hurst', 'threshold', R.FASTER_HurstThreshold};
        if ~isempty(R.FASTER_RefChan), args = [args, {'refchan', R.FASTER_RefChan}]; end
        if ~isempty(R.FASTER_Bandpass), args = [args, {'bandpass', R.FASTER_Bandpass}]; end
        [~, idx] = FASTER_rejchan(EEG, args{:});
        if ~isempty(R.FASTER_RefChan)
            idx = setdiff(idx, R.FASTER_RefChan);
        end
        Bad.Hurst = idx(:)';
        used{end+1} = 'FASTER_Hurst'; doPlot(R.PlotFcn, EEG, Bad.Hurst, R.LogPath, 'Hurst');
    else
        Bad.Hurst = [];
    end

    % 7) CleanRaw: Flatline
    if R.CleanRaw_Flatline
        idx = cleanraw_rejchan(EEG, 'elec', R.IdxDetect, ...
            'measure','flatline', ...
            'threshold', R.Flatline_Sec, ...
            'highpass', R.CleanDrift_Band);
        Bad.Flatline = idx(:)';
        used{end+1} = 'Flatline'; doPlot(R.PlotFcn, EEG, Bad.Flatline, R.LogPath, 'Flatline');
    else
        Bad.Flatline = [];
    end

    % 8) CleanRaw: Noisy channels (corr + line)
    if R.CleanRaw_Noise
        idx = cleanraw_rejchan(EEG, 'elec', R.IdxDetect, ...
            'measure','CleanChan', ...
            'chancorr_crit', R.CleanChan_Corr, ...
            'line_crit',     R.CleanChan_Line, ...
            'max_broken_time', R.CleanChan_MaxBad, ...
            'min_corr_samples', R.CleanChan_NSamp, ...
            'highpass', R.CleanDrift_Band);
        Bad.CleanChan = idx(:)';
        used{end+1} = 'CleanChan'; doPlot(R.PlotFcn, EEG, Bad.CleanChan, R.LogPath, 'CleanChan');
    else
        Bad.CleanChan = [];
    end

    % 9) Known bad (optional)
    Bad.Known = unique(R.KnownBadidx(:)','stable');

    % --------- Combine & summarize ----------
    fields = {'Kurt','Spec','Prob','MeanCorr','Variance','Hurst','Flatline','CleanChan','Known'};
    allBad = [];
    summary = struct();
    for k = 1:numel(fields)
        f = fields{k};
        if ~isfield(Bad,f) || isempty(Bad.(f)), Bad.(f) = []; end
        summary.(f) = numel(Bad.(f));
        allBad = [allBad, Bad.(f)]; %#ok<AGROW>
    end
    % unique while preserving first occurrence order
    Bad.all = unique(allBad, 'stable');
    summary.Total = numel(Bad.all);

    % --------- Action: stash | remove | flag ----------
    switch lower(R.Action)
        case 'remove'
            if ~isempty(Bad.all)
                EEG = pop_select(EEG, 'rmchannel', Bad.all);
                EEG = eeg_checkset(EEG);
            end
        case 'flag'
            % 1=keep, 0=bad; initialize all as 1
            mask = true(1, EEG.nbchan);
            mask(Bad.all(Bad.all>=1 & Bad.all<=EEG.nbchan)) = false;
            EEG.etc.clean_channel_mask = mask(:)';
    end

    % --------- Bookkeeping in EEG.etc ----------
    out = struct();
    out.Bad = Bad;
    out.summary = summary;
    out.detectors_used = used;
    out.IdxDetect = R.IdxDetect;

    EEG.etc.EEGdojo = ensureStructField(EEG.etc, 'EEGdojo', struct());
    EEG.etc.EEGdojo.BadChanIdx = Bad.all;
    EEG.etc.EEGdojo.BadChanLabel = idx2chans(EEG, Bad.all);
    EEG.etc.EEGdojo.BadChanSummary = summary;
    EEG.etc.EEGdojo.BadDetectorsUsed = used;

end

% ---------- local helpers (minimal; no fallbacks) ----------
function doPlot(plotFcn, EEG, idx, logPath, tag)
    if isempty(plotFcn) || isempty(idx), return; end
    try
        plotFcn(EEG, idx, logPath, tag);
    catch
        % silent; plotting should not break pipeline
    end
end

function S = ensureStructField(S, f, v)
    if ~isfield(S, f), S.(f) = v; end
end

