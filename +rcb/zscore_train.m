function [Xz, mu, sd] = zscore_train(X)
% ZSCORE_TRAIN  Z-scores data and returns transformation parameters.
%
%   [Xz, mu, sd] = zscore_train(X)
%
% Inputs:
%   X  - [n x p] data matrix (samples x features).
%
% Outputs:
%   Xz - [n x p] z-scored data.
%   mu - [1 x p] mean of each column.
%   sd - [1 x p] standard deviation of each column.
%
% See also: rcb.zscore_apply

mu = mean(X, 1, 'omitnan');
sd = std(X, 0, 1, 'omitnan');

% Avoid division by zero for constant columns
sd(sd < eps) = 1;

Xz = (X - mu) ./ sd;

end
