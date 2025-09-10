function [S, meta] = behav_similarity(B, varargin)
% BEHAV_SIMILARITY  Subject-by-subject behavioral similarity for IS-RSA.
%
%   [S, meta] = behav_similarity(B, 'Name',Value,...)
%
% PURPOSE
%   Build an N×N similarity matrix from behavioral data for Inter-Subject RSA.
%   - For a single trait per subject (N×1), use rank-based models recommended
%     by the IS-RSA literature: 'nn', 'annak-mean', 'annak-min'.
%   - For multi-item/feature data per subject (N×P), use embedding similarity
%     (correlation or cosine).
%
% INPUT
%   B : [N x 1] vector of a single trait/score per subject
%       or [N x P] matrix of P features (e.g., questionnaire items).
%
% NAME-VALUE OPTIONS
%   'model'   : 'nn' (default) | 'annak-mean' | 'annak-min' | ...
%               'embedding-corr' | 'embedding-cosine'
%               - 'nn'          : nearest-neighbor (similarity decays with rank distance)
%               - 'annak-mean'  : Anna-Karenina (mean of normalized ranks)
%               - 'annak-min'   : Anna-Karenina (min of normalized ranks)
%               - 'embedding-*' : similarity across multi-item vectors
%   'rank'    : true (default) | false
%               For univariate NON-rank models (future extensions) controls ranking.
%               NOTE: For 'nn'/'annak-*' we ALWAYS use ranks (flag is ignored).
%   'zscore'  : true | false (default false)
%               Z-score columns for embedding models before similarity.
%
% OUTPUT
%   S    : [N x N] symmetric similarity matrix (diag = 1)
%   meta : struct with fields:
%          .N, .P, .model, .rank, .zscore
%
% DESIGN NOTES (IS-RSA best practices)
%   - Univariate trait → rank transform, then NN or Anna-K models.
%   - Keep matrices UNSORTED for stats; sort only for visualization.
%   - Use Spearman and Mantel-style permutations downstream for IS-RSA.
%
% EXAMPLES
%   % 1) Univariate trait (NN model)
%   S = behav_similarity(trait, 'model','nn');
%
%   % 2) Anna-K mean-rank model
%   S = behav_similarity(trait, 'model','annak-mean');
%
%   % 3) Multi-item embedding with cosine similarity (and z-scoring)
%   S = behav_similarity(items, 'model','embedding-cosine', 'zscore', true);
%
% v1.1 (2025) – aligned to IS-RSA tutorial conventions; corrected rank handling.

% ------------------ Parse & validate ------------------
if ndims(B) ~= 2
    error('behav_similarity:InputDim', 'B must be 2-D [N x P].');
end
[N, P] = size(B);

ip = inputParser;
ip.addParameter('model','nn');
ip.addParameter('rank',true, @(x)islogical(x)&&isscalar(x));
ip.addParameter('zscore',false, @(x)islogical(x)&&isscalar(x));
ip.parse(varargin{:});
opt   = ip.Results;
model = lower(string(opt.model));

% ------------------ Univariate branch (P == 1) ------------------
if P == 1
    x = B(:);
    if any(isnan(x))
        error('behav_similarity:NaN', ...
              'NaNs found in univariate behavioral vector. Impute/remove before calling.');
    end

    % For NN / Anna-K, ALWAYS rank (ignore user flag); warn if user set rank=false.
    isRankModel = any(strcmpi(model, {"nn","annak-mean","annak-min"}));
    if isRankModel && ~opt.rank
        warning('behav_similarity:RankIgnored', ...
            ['For models ''nn''/''annak-*'', rank-transform is required and will be applied ', ...
             '(the ''rank'' option is ignored).']);
    end

    % Decide representation r (ranked or raw)
    if isRankModel
        r = tiedrank(x);                  % 1..N (ties averaged)
    else
        if opt.rank
            r = tiedrank(x);
        else
            r = x(:);
        end
    end

    % Map ranks to [0,1] when needed
    if isRankModel
        if N == 1
            r01 = 1;  % degenerate case
        else
            r01 = (r - 1) / (N - 1);      % 0..1
        end
    end

    switch model
        case "nn"
            % Nearest-neighbor: similarity decays with rank distance
            if N == 1
                S = 1;
            else
                S = 1 - abs(r - r.') / (N - 1);
                S(1:N+1:end) = 1;
            end

        case "annak-mean"
            % Similarity increases with absolute (normalized) rank
            S = (r01 + r01.') / 2;
            S(1:N+1:end) = 1;

        case "annak-min"
            % Stricter variant: min of normalized ranks
            S = min(r01, r01.');
            S(1:N+1:end) = 1;

        otherwise
            error('behav_similarity:model', ...
                'Model "%s" requires univariate input and must be nn/annak-mean/annak-min.', model);
    end

% ------------------ Embedding branch (P >= 2) ------------------
else
    X = B;

    if any(isnan(X(:)))
        error('behav_similarity:NaN', ...
              'NaNs found in behavioral matrix. Impute/remove before calling.');
    end

    if opt.zscore
        % Z-score columns; guard zero-variance columns
        mu = mean(X,1);
        sd = std(X,[],1);
        sd(sd==0) = 1;  % avoid division by zero; constant columns become zeros after centering
        X = (X - mu) ./ sd;
    end

    switch model
        case "embedding-corr"
            % subjects×subjects correlation
            S = corrcoef(X.');  % corrcoef expects variables in columns
            % corrcoef sets diagonal to 1

        case "embedding-cosine"
            % Cosine similarity across subject vectors
            d = sqrt(sum(X.^2,2)); d(d==0) = eps;
            Xn = X ./ d;
            S = Xn * Xn.';                 % symmetric, diag ~ 1
            S(1:N+1:end) = 1;

        otherwise
            error('behav_similarity:model', ...
                 'Unknown embedding model "%s". Use embedding-corr or embedding-cosine.', model);
    end
end

% ------------------ Output meta ------------------
meta = struct('N',N,'P',P,'model',char(model), ...
              'rank',logical(opt.rank), 'zscore',logical(opt.zscore));
end
