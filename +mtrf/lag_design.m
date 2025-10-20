function [Xlag, kept_idx, lag_samp] = lag_design(X, lags, varargin)
% LAG_DESIGN  Build a time-lagged design matrix (mTRF-style).
%
%   [Xlag, kept_idx, lag_samp] = lag_design(X, lags, 'Name', Value, ...)
%
% Inputs
%   X      : [T x F] feature matrix (rows=time/samples; cols=features).
%   lags   : vector of lags (ms or samples; see 'units').
%            Convention (mTRF): positive lag uses X(t - lag) to predict Y(t).
%
% Name-Value options (defaults in brackets)
%   'units'   : 'samples' | 'ms'                                  ['samples']
%   'fs'      : sampling rate (Hz) if 'units'=='ms'               []
%   'zeropad' : logical: keep T rows and pad edges, else trim      [true]
%   'pad'     : 'zero' | 'nan' | 'edge'  (only if zeropad=true)   ['zero']
%               - 'zero' -> pad with zeros
%               - 'nan'  -> pad with NaNs (drop later before CV)
%               - 'edge' -> replicate the nearest valid row
%   'bias'    : logical: prepend a column of ones                   [false]
%
% Outputs
%   Xlag     : [T x (F*L)] if zeropad==true, else [T_keep x (F*L)]
%              where L = numel(lags). Column blocks are ordered by lags(:).
%   kept_idx : indices (in original X) of rows kept in Xlag.
%              If zeropad==true -> 1:T; if zeropad==false -> trimmed range.
%   lag_samp : integer lags in samples (same order as input lags).
%
% Notes
%   • With zeropad==false, rows that would require out-of-bounds samples
%     are removed (contiguous block). This matches mTRF’s no-padding mode.
%   • With zeropad==true and 'pad'='nan', remember to drop rows with NaNs
%     before fitting or ensure your CV splitter masks them.
%   • Column order: [X at lags(1), X at lags(2), ...], each block has F cols.
%
% Example (ms lags, zero-pad with edge replication, include bias):
%   [Xlag, idx, lS] = lag_design(X, -100:10:500, 'units','ms','fs',64, ...
%                                'zeropad',true,'pad','edge','bias',true);
%
% Example (sample lags, trim):
%   [Xlag, idx, lS] = lag_design(X, -8:8, 'units','samples', ...
%                                'zeropad',false,'bias',false);

ip = inputParser;
ip.addParameter('units', 'samples', @(s) any(strcmpi(s,{'samples','ms'})));
ip.addParameter('fs', [], @(x) isempty(x) || (isscalar(x) && x>0));
ip.addParameter('zeropad', true, @islogical);
ip.addParameter('pad', 'zero', @(s) any(strcmpi(s,{'zero','nan','edge'})));
ip.addParameter('bias', false, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

[T, F] = size(X);
lags = lags(:);
L = numel(lags);

% Convert lags to integer samples
switch lower(opt.units)
    case 'samples'
        lag_samp = round(lags);
    case 'ms'
        assert(~isempty(opt.fs), 'lag_design: ''fs'' is required when units="ms".');
        lag_samp = round(lags .* opt.fs ./ 1000);
    otherwise
        error('lag_design: unknown units "%s".', opt.units);
end

% Determine trimmed range if not zero-padding
if ~opt.zeropad
    maxPos = max(0, max(lag_samp));      % largest positive lag
    maxNeg = max(0, -min(lag_samp));     % absolute of most negative lag
    kept_idx = (1+maxPos) : (T - maxNeg);
    if isempty(kept_idx)
        error('lag_design: lags range exceeds available samples.');
    end
    Tkeep = numel(kept_idx);
    Xlag = zeros(Tkeep, F*L, 'like', X);
else
    kept_idx = (1:T).';
    Xlag = zeros(T, F*L, 'like', X);     % will overwrite below
end

% Precompute padding factories
use_edge = strcmpi(opt.pad,'edge');
use_nan  = strcmpi(opt.pad,'nan');

% Build each lag block
for i = 1:L
    s = lag_samp(i);
    cols = ( (i-1)*F + 1 ) : (i*F);

    if s >= 0
        % Xlag(t,:) uses X(t - s,:)  -> pad at TOP when zeropad
        core = X(1:T-s, :);
        if opt.zeropad
            if s > 0
                if use_edge
                    topPad = repmat(X(1,:), s, 1);
                elseif use_nan
                    topPad = nan(s, F, class(X));
                else
                    topPad = zeros(s, F, class(X));
                end
                Xi = [topPad; core];
            else
                Xi = core;
            end
            Xlag(:, cols) = Xi;
        else
            Xi = core( (1+max(0,0)) : end, : ); %#ok<NASGU> (clarity)
            % When trimming, we just take the rows that align with kept_idx:
            Xlag(:, cols) = X(kept_idx - s, :);
        end

    else
        % Negative lag: s<0 => uses future X(t + |s|); pad at BOTTOM
        s2 = -s;
        core = X(1+s2:T, :);
        if opt.zeropad
            if s2 > 0
                if use_edge
                    botPad = repmat(X(end,:), s2, 1);
                elseif use_nan
                    botPad = nan(s2, F, class(X));
                else
                    botPad = zeros(s2, F, class(X));
                end
                Xi = [core; botPad];
            else
                Xi = core;
            end
            Xlag(:, cols) = Xi;
        else
            Xlag(:, cols) = X(kept_idx - s, :);  % since s is negative
        end
    end
end

% Add bias if requested
if opt.bias
    if opt.zeropad
        Xlag = [ones(size(Xlag,1),1, 'like', X), Xlag];
    else
        Xlag = [ones(numel(kept_idx),1, 'like', X), Xlag];
    end
end
end
