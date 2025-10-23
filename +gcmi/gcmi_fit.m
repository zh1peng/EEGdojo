function out = gcmi_fit(Y, X, Cov, varargin)
% GCMI_FIT Temporal (conditional) mutual information across lags with null inference.
%
% OUT = GCMI_FIT(Y, X, Cov, 'Name',Value,...)
% OUT = GCMI_FIT(Y, X) % no covariates (plain MI)
%
% Inputs (rows = time samples)
% Y : [T x Dy] neural time series (Dy=1 typical; multivariate OK)
% X : [T x Dx] stimulus feature(s) of interest (one or many)
% Cov : [T x Dz] OPTIONAL continuous covariates to condition on
% (e.g., luminance, motion, audio); pass [] to omit.
%
% Name-Value options (sane defaults in brackets)
% 'Lags' : vector of lags (ms or samples, see 'Units') [ -500:10:1000 ]
% 'Units' : 'ms' | 'samples' [ 'ms' ]
% 'Fs' : sampling rate in Hz, required if Units='ms' or NullMethod='block' [ [] ]
% 'NullMethod' : 'circular' | 'block' | 'phase' [ 'circular' ]
% 'Nperm' : number of permutations for null [ 1000 ]
% 'BlockLenSec' : block length (s) for 'block' null [ 10 ]
% 'BiasCorrect' : true/false, Gaussian-entropy bias correction [ false ]
% 'RngSeed' : scalar RNG seed [ 1 ]
% 'KeepNull' : store full null MI matrix in OUT.null_MI [ true ]
% 'Verbose' : print progress [ true ]
%
% Outputs (all per-lag vectors are column-shaped)
% out.MI : [nLag x 1] MI or CMI (bits) at each lag
% out.z : [nLag x 1] z-score vs. permutation null
% out.p : [nLag x 1] upper-tail p-value vs. null (maximization)
% out.null_mean : [nLag x 1] null mean (bits)
% out.null_std : [nLag x 1] null std (bits)
% out.null_MI : [nPerm x nLag] null MI per lag (if KeepNull=true)
% out.lags : [nLag x 1] lags in input units
% out.lag_samples : [nLag x 1] lags in samples
% out.lags_ms : [nLag x 1] lags in ms (if Fs provided)
% out.valid_nsamps : [nLag x 1] number of overlapping samples used
% out.info : struct of settings & shapes (Units, Fs, NullMethod, etc.)
%
% Requirements
% - GCMI toolbox on MATLAB path (robince/gcmi): copnorm, mi_gg, cmi_ggg, etc.
% - Inputs must have the same number of rows (time samples).
%
% Notes
% - We trim to overlapping samples for each lag (no padding).
% - For nulls we misalign/shuffle X ONLY (variable of interest) and keep Y and Cov fixed.
% - The CMI calculated is I(X(t-lag); Y(t) | Cov(t-lag)).
% - Use 'phase' null to preserve spectrum; 'block' to respect nonstationarity; 'circular' as fast default.
%
% Example
% out = gcmi_fit(Y, X, Cov, 'Lags',-400:10:800,'Units','ms','Fs',128, ...
% 'NullMethod','circular','Nperm',1000);
%
% % Unique info in X about Y beyond Cov is out.MI (bits) per lag, with z and p.

% ---------------------- parse & validate ----------------------
if nargin < 3, Cov = []; end
ip = inputParser;
ip.addParameter('Lags', -200:10:500);
ip.addParameter('Units', 'ms', @(s) any(strcmpi(s,{'ms','samples'})));
ip.addParameter('Fs', []);
ip.addParameter('NullMethod', 'circular', @(s) any(strcmpi(s,{'circular','block','phase'})));
ip.addParameter('Nperm', 1000, @(x) isscalar(x) && x>=0);
ip.addParameter('BlockLenSec', 10, @(x) isscalar(x) && x>0);
ip.addParameter('BiasCorrect', false, @islogical);
ip.addParameter('RngSeed', 1, @(x) isscalar(x) && isnumeric(x));
ip.addParameter('KeepNull', true, @islogical);
ip.addParameter('Verbose', true, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

% --- Validate dimensions ---
[T, Dy] = size(Y);
[Tx, Dx] = size(X);
if Tx ~= T, error('Y and X must have the same number of rows (time samples).'); end
hasCov = ~isempty(Cov);
if hasCov
    [Tc, Dz] = size(Cov);
    if Tc ~= T, error('Y, X, and Cov must share the same number of rows.'); end
else
    Dz = 0;
end

% --- Lags in samples & ms ---
switch lower(opt.Units)
case 'samples'
    lag_samples = round(opt.Lags(:));
    if ~isempty(opt.Fs)
        lags_ms = (lag_samples(:) / opt.Fs) * 1000;
    else
        lags_ms = nan(numel(lag_samples),1);
    end
case 'ms'
    if isempty(opt.Fs), error('Fs is required when Units="ms".'); end
    lag_samples = round(opt.Lags(:) * opt.Fs / 1000);
    lags_ms = opt.Lags(:);
otherwise
    error('Unknown Units: %s', opt.Units);
end
nLag = numel(lag_samples);

% --- Minimum samples check ---
% GCMI requires N > D (total dimensions) for cov matrix to be non-singular
if hasCov
    minSamps = Dx + Dy + Dz;
else
    minSamps = Dx + Dy;
end

% ---------------------- prep transforms ----------------------
rng(opt.RngSeed);

% Copula-normalize once (columns independent)
Yc = copnorm(Y);
Xc = copnorm(X);
if hasCov, Zc = copnorm(Cov); else, Zc = []; end

% Precompute per-lag aligned indices (overlapping samples)
[idxY, idxX, idxZ, nOverlap] = build_lag_indices(T, lag_samples, hasCov);

% Allocate outputs
MI = nan(nLag,1);
null_mu = nan(nLag,1);
null_sd = nan(nLag,1);
validN = nOverlap(:);
store_null = opt.KeepNull && (opt.Nperm > 0);
if store_null
    nullMI = nan(opt.Nperm, nLag);
else
    nullMI = [];
end

% ---------------------- observed MI/CMI per lag ----------------------
if opt.Verbose, fprintf('[gcmi_fit] Calculating observed MI/CMI...\n'); end
for i = 1:nLag
    % PATCH: Skip if N <= D (not just N <= 1)
    if validN(i) <= minSamps, continue; end
    
    y = Yc(idxY{i}, :);
    x = Xc(idxX{i}, :);
    if hasCov
        z = Zc(idxZ{i}, :);
        MI(i) = cmi_ggg(x, y, z, opt.BiasCorrect); % bits
    else
        MI(i) = mi_gg(x, y, opt.BiasCorrect); % bits
    end
end

% ---------------------- permutation null ----------------------
if opt.Nperm > 0
    switch lower(opt.NullMethod)
    case 'circular'
        % choose offsets that are at least 10% of T away from 0 shift
        minGap = max(1, round(0.10*T));
        validOff = (minGap : (T - minGap));
        if isempty(validOff) % Handle very short T
             validOff = 1:max(1, T-1);
        end
    case 'block'
        % PATCH: 'block' method requires Fs to interpret BlockLenSec
        if isempty(opt.Fs)
            error('Fs is required for NullMethod="block" (to interpret BlockLenSec).');
        end
        blockLenSamples = round(opt.BlockLenSec * opt.Fs);
        if blockLenSamples < 1, blockLenSamples = 1; end
        if blockLenSamples > T, blockLenSamples = T; end
        blk = make_blocks(T, blockLenSamples);
    case 'phase'
        % precompute phase-scrambled versions of Xc per permutation (col-wise)
        % (done inside loop for low memory; deterministic with RngSeed)
    otherwise
        error('Unknown NullMethod.');
    end

    for p = 1:opt.Nperm
        % Build permuted Xc only (misalign variable of interest)
        switch lower(opt.NullMethod)
            case 'circular'
                off = validOff(randi(numel(validOff)));
                Xp = circshift(Xc, off, 1);
            case 'block'
                order = randperm(numel(blk));
                % PATCH: Inlined blk_apply logic
                Xp_parts = cellfun(@(ix) Xc(ix,:), blk(order), 'uni', false);
                Xp = vertcat(Xp_parts{:});
            case 'phase'
                Xp = phase_randomize_cols(Xc);
        end

        % MI per lag using permuted X
        for i = 1:nLag
            % PATCH: Skip if N <= D
            if validN(i) <= minSamps, continue; end
            
            y = Yc(idxY{i}, :);
            x = Xp(idxX{i}, :); % Use permuted X, but same lag indices
            if hasCov
                z = Zc(idxZ{i}, :); % Covariates kept fixed (aligned to Y/Xp)
                val = cmi_ggg(x, y, z, opt.BiasCorrect);
            else
                val = mi_gg(x, y, opt.BiasCorrect);
            end
            if store_null, nullMI(p,i) = val; end
        end

        if opt.Verbose && mod(p, max(1,floor(opt.Nperm/10)))==0
            fprintf('[gcmi_fit] Permutations: %d / %d\n', p, opt.Nperm);
        end
    end

    % Null moments & z/p per lag (upper tail)
    null_mu = mean(nullMI, 1, 'omitnan')';
    null_sd = std( nullMI,  [], 1, 'omitnan')';

else
    null_mu(:) = NaN; null_sd(:) = NaN;
end

% z and p (upper-tail, permutation-based)
z = nan(nLag,1);
p = nan(nLag,1);
if opt.Nperm > 0
    z = (MI - null_mu) ./ (null_sd + eps);
    % permutation p: proportion of null >= observed (with +1 rule)
    for i = 1:nLag
        % PATCH: Skip if N <= D
        if validN(i) <= minSamps, p(i)=NaN; continue; end
        if store_null
            p(i) = (1 + sum(nullMI(:,i) >= MI(i))) / (opt.Nperm + 1);
        else
            p(i) = NaN; % null not stored; user has z, null mean/std
        end
    end
end

% ---------------------- pack output ----------------------
out.MI = MI;
out.z = z;
out.p = p;
out.null_mean = null_mu;
out.null_std = null_sd;
out.null_MI = nullMI;
out.lags = opt.Lags(:);
out.lag_samples = lag_samples(:);
out.lags_ms = lags_ms(:);
out.valid_nsamps = validN(:);
out.info = struct( ...
    'Units', opt.Units, ...
    'Fs', opt.Fs, ...
    'NullMethod', opt.NullMethod, ...
    'Nperm', opt.Nperm, ...
    'BlockLenSec', opt.BlockLenSec, ...
    'BiasCorrect', opt.BiasCorrect, ...
    'RngSeed', opt.RngSeed, ...
    'Dy', Dy, ...
    'Dx', Dx, ... % PATCH: Added Dx
    'Dz', Dz ...  % PATCH: Added Dz
);

if opt.Verbose
    fprintf('[gcmi_fit] Done. %d lags; Nperm=%d; Null=%s\n', nLag, opt.Nperm, opt.NullMethod);
end

end % -- gcmi_fit --

% ===================== helper functions =====================

function [idxY, idxX, idxZ, nOverlap] = build_lag_indices(T, lag_samples, hasCov)
% Build per-lag overlapping indices for Y, X, (Cov), trimming edges.
% Note: X and Z (Cov) are aligned: I(X(t-s); Y(t) | Z(t-s))
nLag = numel(lag_samples);
idxY = cell(nLag,1);
idxX = cell(nLag,1);
if hasCov, idxZ = cell(nLag,1); else, idxZ = []; end
nOverlap = zeros(nLag,1);

for i = 1:nLag
    s = lag_samples(i);
    if s >= 0
        % s = +10. Pair Y(11:T) with X(1:T-10)
        % Y(t) vs X(t-s)
        iy = (1+s):T;
        ix = 1:(T-s);
    else
        % s = -10. Pair Y(1:T-10) with X(11:T)
        % Y(t) vs X(t-s)
        s2 = -s;
        iy = 1:(T-s2);
        ix = (1+s2):T;
    end
    idxY{i} = iy(:);
    idxX{i} = ix(:);
    if hasCov, idxZ{i} = ix(:); end % Z is aligned with X
    nOverlap(i) = numel(iy);
end
end

function blk = make_blocks(T, blockLen)
% Return cell array of index vectors for contiguous blocks (last may be shorter).
if blockLen < 1, blockLen = 1; end
nblk = ceil(T / blockLen);
blk = cell(1, nblk);
for b = 1:nblk
    a = (b-1)*blockLen + 1;
    z = min(b*blockLen, T);
    blk{b} = a:z;
end
end

function Xp = phase_randomize_cols(X)
% Column-wise phase randomization preserving power spectrum; returns real signal.
[T, D] = size(X);
Xp = zeros(T, D, 'like', X);

for j = 1:D
    x  = X(:,j);
    xf = fft(x);
    N  = numel(xf);
    xf_new = xf;

    % PATCH: Correctly identify positive frequencies to randomize
    % (exclude DC and Nyquist)
    nPos = floor((N-1)/2); % Number of positive freqs
    pos = 2:(nPos + 1);    % Indices of positive freqs
    
    if isempty(pos), continue; end % Skip if T is too small

    % Random phases
    ph = exp(1i * 2*pi*rand(numel(pos),1));

    % Apply to positive freqs, preserve magnitudes
    xf_new(pos) = abs(xf(pos)) .* ph;

    % Enforce conjugate symmetry for a real time-domain signal
    neg = N - pos + 2; % Indices of negative freqs (in reverse order)
    xf_new(neg) = conj(xf_new(pos)); % This is correct, no flip needed

    % Keep DC (1) and Nyquist (if even length: N/2+1) unchanged
    Xp(:,j) = real(ifft(xf_new));
end

% Standardize (optional; copula ranks unaffected, but keeps scale tidy)
% Xp = (Xp - mean(Xp,1)) ./ (std(Xp,0,1) + eps);
end