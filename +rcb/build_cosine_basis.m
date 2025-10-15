function B = build_cosine_basis(lags_ms, centers_ms, width_ms)
% BUILD_COSINE_BASIS Generates a basis set of overlapping raised cosines.
%
%   B = build_cosine_basis(lags_ms, centers_ms, width_ms)
%
%   This function creates a set of raised cosine functions that span a specified
%   temporal domain (lags). These are often used as a basis set to smoothly
%   model a temporal response function (TRF) with fewer parameters than a
%   standard finite impulse response (FIR) model.
%
%   Each basis function is a Hanning window, defined as:
%   y(t) = 0.5 * (1 + cos(pi * (t - c) / w)) for |t - c| <= w
%
% Inputs:
%   lags_ms    - [1 x L] or [L x 1] vector of time points for the lags.
%   centers_ms - [1 x K] or [K x 1] vector specifying the center of each cosine basis.
%   width_ms   - Scalar defining the full width of each cosine bell.
%
% Output:
%   B          - [L x K] matrix where each column is a basis function, L2-normalized.
%
% See also: rcb.plot_cosine_basis

% Input validation
assert(isvector(lags_ms), 'lags_ms must be a vector.');
assert(isvector(centers_ms), 'centers_ms must be a vector.');
assert(isscalar(width_ms) && width_ms > 0, 'width_ms must be a positive scalar.');

t = lags_ms(:);         % Ensure t is a column vector [L x 1]
C = centers_ms(:)';     % Ensure C is a row vector [1 x K]
w = width_ms;

% Calculate the phase of each lag relative to each center
% phi will be an [L x K] matrix
phi = (t - C) ./ w;

% Create the raised cosine bells
% This is equivalent to a Hanning window shape
B = 0.5 * (1 + cos(pi * max(-1, min(1, phi))));

% Set values outside the window width to zero
B(abs(phi) > 1) = 0;

% Normalize each basis function (column) to have a unit L2-norm
B = B ./ max(eps, sqrt(sum(B.^2, 1)));

end