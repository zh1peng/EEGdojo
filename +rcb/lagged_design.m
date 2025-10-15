function [Xlag, lag_samp] = lagged_design(X, lags, varargin)
% RCB.LAGGED_DESIGN  Build a classic lagged design matrix from features X.
%
%   [Xlag, lag_samp] = rcb.lagged_design(X, lags, 'fs', fs, 'units','ms')
%   [Xlag, lag_samp] = rcb.lagged_design(X, lags, 'units','samples')
%
% Inputs
%   X      : [T x F] feature matrix.
%   lags   : vector of lags (ms or samples, see 'units').
%            Convention: positive lag = past predictors -> present y.
%            I.e., a lag of +L samples uses X(t-L) to predict y(t).
%
% Name-Value Pairs (optional)
%   'fs'     : sampling rate (Hz). Required if 'units'='ms'. (default [])
%   'units'  : 'ms' (default) or 'samples'
%   'pad'    : edge padding for impossible shifts: 'nan' (default), 'zero', or 'edge'
%
% Outputs
%   Xlag     : [T x (F * numel(lags))] concatenated lagged features.
%              Column blocks are ordered by lags(:) in that order.
%   lag_samp : [numel(lags) x 1] integer lags in samples (rounded).
%
% Notes
%   • Rows where a shift is impossible are padded per 'pad'.
%   • If you keep 'pad'='nan', remember to drop rows with any NaNs before CV.
%   • Pair this with rcb_fit(...,'DesignType','lagged') to get TRFs by reshaping
%     weights to [L x F].

ip = inputParser;
ip.addParameter('fs', [], @(x) isempty(x) || (isscalar(x) && x>0));
ip.addParameter('units', 'ms', @(s) any(strcmpi(s,{'ms','samples'})));
ip.addParameter('pad', 'nan', @(s) any(strcmpi(s,{'nan','zero','edge'})));
ip.parse(varargin{:});
opt = ip.Results;

[T, F] = size(X);
lags = lags(:);                   % column vector
Lnum  = numel(lags);

% Convert to samples
switch lower(opt.units)
    case 'ms'
        assert(~isempty(opt.fs), 'rcb.lagged_design: fs is required when units="ms".');
        lag_samp = round(lags .* opt.fs ./ 1000);
    case 'samples'
        lag_samp = round(lags);   % ensure integer shifts
    otherwise
        error('Unknown units: %s', opt.units);
end

% Preallocate and select padding behavior
switch lower(opt.pad)
    case 'nan'
        padfun = @(n) nan(n, F, class(X));
    case 'zero'
        padfun = @(n) zeros(n, F, class(X));
    case 'edge'
        % Will replicate the first/last valid row as needed
        padfun = []; % handled inline
    otherwise
        error('Unknown pad: %s', opt.pad);
end

Xlag = nan(T, F*Lnum, class(X));    % initialize

for i = 1:Lnum
    s = lag_samp(i);
    colIdx = ( (i-1)*F + 1 ) : (i*F);

    if s >= 0
        % Use X(t-s) to predict y(t): NaNs (or pad) at the TOP
        core = X(1:T-s, :);
        if strcmpi(opt.pad,'edge')
            topPad = repmat(X(1,:), s, 1);
        else
            topPad = padfun(s);
        end
        Xi = [topPad; core];
    else
        % Negative lag: uses future predictors X(t+|s|); pad at the BOTTOM
        s2 = -s;
        core = X(1+s2:T, :);
        if strcmpi(opt.pad,'edge')
            botPad = repmat(X(end,:), s2, 1);
        else
            botPad = padfun(s2);
        end
        Xi = [core; botPad];
    end

    Xlag(:, colIdx) = Xi;
end
end
