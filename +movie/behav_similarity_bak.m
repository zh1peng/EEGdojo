function [S, meta] = behav_similarity(B, varargin)
% BEHAV_SIMILARITY  Subject-by-subject behavioral similarity for IS-RSA.
%
%   [S, meta] = behav_similarity(B, 'Name',Value,...)
%
% PURPOSE
%   Build an N×N similarity matrix from behavioral data for Inter-Subject RSA.
%   - Univariate (N×1): rank-based models: 
%       'nn'|'nearest-neighbor'|'nn-absdiff'  : 1 - |r_i - r_j|/(N-1)
%       'divergence'                          : 1 - mean(r01_i, r01_j)
%       'convergence' or 'annak-min'          : min(r01_i, r01_j)
%       'annak-mean'                          : mean(r01_i, r01_j)
%   - Multivariate (N×P, P>1): embedding similarity 
%       'embedding-corr' | 'embedding-cosine'
%
% NAME-VALUE OPTIONS
%   'model'   : see above. Default 'nn'.
%   'rank'    : true|false (ignored for univariate models; ranks always used there)
%   'zscore'  : true|false (default false; for embedding models)
%
% OUTPUT
%   S    : [N x N] symmetric similarity (diag = 1), in [0,1] for univariate models
%   meta : struct with fields: .N, .P, .model, .rank, .zscore
%
% DESIGN NOTES
%   - Univariate: rank-transform (tiedrank), then compute pairwise similarity.
%   - Keep S unsorted for stats; only sort for visualization.
%   - Downstream IS-RSA: use Spearman/Mantel-style permutations on upper triangles.
% 
% EXAMPLES
% Univariate maturity score (higher = older)
% S_nn   = behav_similarity(maturity, 'model','nn');            % nearest-neighbor
% S_div  = behav_similarity(maturity, 'model','divergence');    % younger pairs similar
% S_conv = behav_similarity(maturity, 'model','convergence');   % older pairs similar

% Multivariate questionnaire (rows=subjects, cols=items)
% S_cos  = behav_similarity(items, 'model','embedding-cosine','zscore',true);

% v1.3 (2025) – add 'divergence' and 'convergence' (paper-aligned aliases); 
%               normalize univariate similarities to [0,1]; small refactors.

% ------------------ Parse & validate ------------------
if ndims(B) ~= 2, error('behav_similarity:InputDim','B must be 2-D [N x P].'); end
[N, P] = size(B);

ip = inputParser;
ip.addParameter('model','nn');
ip.addParameter('rank',true, @(x)islogical(x)&&isscalar(x));
ip.addParameter('zscore',false, @(x)islogical(x)&&isscalar(x));
ip.parse(varargin{:});
opt = ip.Results;

% Normalize model name
model = lower(string(opt.model));
model = replace(model, "_", "-");
model_char = char(model);

% ------------------ Univariate branch (P == 1) ------------------
if P == 1
    x = B(:);
    if any(isnan(x))
        error('behav_similarity:NaN', ...
              'NaNs in behavioral vector. Impute/remove before calling.');
    end

    % For univariate models, we ALWAYS use ranks
    r  = tiedrank(x);                  % 1..N (ties averaged)
    if N == 1
        r01 = 1;                       % normalized ranks 0..1
    else
        r01 = (r - 1) / (N - 1);
    end

    switch model_char
        case {'nn','nearest-neighbor','nn-absdiff'}
            % Nearest-neighbor: similarity decays with rank distance
            if N == 1
                S = 1;
            else
                S = 1 - abs(r - r.') / (N - 1);  % in [0,1]
                S(1:N+1:end) = 1;
            end

        case {'divergence'}
            % Younger pairs more similar: 1 - mean(r01_i, r01_j)
            % Paper: "sample maximum minus pair average"; with max=1 → 1 - mean(...)
            S = 1 - (r01 + r01.')/2;
            S(1:N+1:end) = 1;

        case {'convergence','annak-min'}
            % Older pairs more similar: pair minimum of normalized ranks
            S = min(r01, r01.');
            S(1:N+1:end) = 1;

        case {'annak-mean'}
            % Mean of normalized ranks (older-heavier average)
            S = (r01 + r01.')/2;
            S(1:N+1:end) = 1;

        otherwise
            error('behav_similarity:model', ...
                'Univariate model "%s" must be one of: nn, divergence, convergence, annak-mean, annak-min.', ...
                model_char);
    end

% ------------------ Multivariate branch (P >= 2) ------------------
else
    X = B;
    if any(isnan(X(:)))
        error('behav_similarity:NaN', ...
              'NaNs in behavioral matrix. Impute/remove before calling.');
    end
    if opt.zscore
        mu = mean(X,1);
        sd = std(X,[],1);
        sd(sd==0) = 1;
        X = (X - mu) ./ sd;
    end

    switch model_char
        case 'embedding-corr'
            S = corrcoef(X.');                 % N×N; diag = 1
        case 'embedding-cosine'
            d = sqrt(sum(X.^2,2)); d(d==0) = eps;
            Xn = X ./ d;
            S = Xn * Xn.';                     % ~[-1,1], diag ~ 1 (exact 1 if no zero rows)
            S(1:N+1:end) = 1;
        otherwise
            error('behav_similarity:model', ...
                 'Unknown embedding model "%s". Use embedding-corr or embedding-cosine.', model_char);
    end
end

% ------------------ Output meta ------------------
meta = struct('N',N,'P',P,'model',model_char, ...
              'rank',true, 'zscore',logical(opt.zscore));
end
