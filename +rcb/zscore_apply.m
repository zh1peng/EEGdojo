function Xz = zscore_apply(X, mu, sd)
% ZSCORE_APPLY  Applies a pre-computed z-score transformation.
%
%   Xz = zscore_apply(X, mu, sd)
%
% Inputs:
%   X  - [m x p] data matrix to transform.
%   mu - [1 x p] mean vector from the training set.
%   sd - [1 x p] standard deviation vector from the training set.
%
% Output:
%   Xz - [m x p] transformed data.
%
% See also: rcb.zscore_train

% Ensure sd vector is the same size as mu
if size(sd, 2) ~= size(mu, 2)
    error('mu and sd must have the same number of columns.');
end

% Avoid division by zero for constant columns
sd(sd < eps) = 1;

Xz = (X - mu) ./ sd;

end
