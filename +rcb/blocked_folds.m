function folds = blocked_folds(T, K, purge_pts)
% BLOCKED_FOLDS  Creates indices for K-fold blocked cross-validation.
%
%   folds = blocked_folds(T, K, purge_pts)
%
%   Partitions a dataset of T samples into K contiguous blocks for
%   cross-validation. A purge buffer can be specified to exclude samples
%   around the boundaries of the test set from the training set, preventing
%   information leakage from temporal dependencies.
%
% Inputs:
%   T         - Total number of samples.
%   K         - Number of folds.
%   purge_pts - Number of samples to purge from the training set on both
%               sides of the test set block.
%
% Output:
%   folds     - [K x 1] struct array with fields:
%               .trainIdx - Logical index for training samples.
%               .testIdx  - Logical index for testing samples.

if nargin < 3, purge_pts = 0; end

assert(K > 1, 'Number of folds (K) must be at least 2.');
assert(T / K >= 1, 'Not enough samples for K folds.');

edges = round(linspace(0, T, K + 1));
folds = repmat(struct('trainIdx', [], 'testIdx', []), K, 1);

for k = 1:K
    t1 = edges(k) + 1;
    t2 = edges(k + 1);
    
    test_mask = false(T, 1);
    test_mask(t1:t2) = true;
    
    train_mask = ~test_mask;
    
    % Purge samples around the test block
    p1 = max(1, t1 - purge_pts);
    p2 = min(T, t2 + purge_pts);
    train_mask(p1:p2) = false;
    
    folds(k).trainIdx = train_mask;
    folds(k).testIdx = test_mask;
end

end
