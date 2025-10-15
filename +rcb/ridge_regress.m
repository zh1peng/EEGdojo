function W = ridge_regress(X, y, lambda, varargin)
% RIDGE_REGRESS  Ridge regression with robust solves and multi-target support.
%
%   W = ridge_regress(X, y, lambda)
%   W = ridge_regress(X, y, lambda, 'reg_mask', ones(size(X,2),1))
%
% Inputs
%   X        : [n x p] predictors
%   y        : [n x m] response(s); m can be 1 or >1
%   lambda   : non-negative scalar
%
% Name-Value
%   'reg_mask' : [p x 1] 0/1 (or weights) for per-feature L2 penalty.
%                Default all ones (penalize everything). Use 0 to skip bias.
%
% Output
%   W : [p x m] weights
%
% Notes
%   • Uses Woodbury for p>n; Cholesky in both branches; SVD fallback.
%   • If you pass a bias column in X and want it unpenalized, set reg_mask(end)=0.

ip = inputParser;
ip.addParameter('reg_mask', [], @(v) isempty(v) || (isnumeric(v) && isvector(v)));
ip.parse(varargin{:});
reg_mask = ip.Results.reg_mask;

assert(isscalar(lambda) && isreal(lambda) && lambda >= 0, 'Lambda must be a non-negative scalar.');

[n, p] = size(X);
if size(y,1) ~= n, error('X and y must have the same number of rows.'); end
if isempty(reg_mask), reg_mask = ones(p,1, 'like', X); end
if numel(reg_mask) ~= p, error('reg_mask must be length p.'); end
reg_mask = reshape(reg_mask, [], 1);

% Make sure y is 2D
if isvector(y), y = reshape(y, n, []); end

% Fast path: p > n -> Woodbury
if p > n
    % K = X*X' + lambda*I  (note: lambda effectively scaled uniformly here)
    % If you provide a non-uniform reg_mask, Woodbury becomes messy;
    % we fallback to the p×p formulation in that case.
    if all(abs(reg_mask - 1) < eps)
        I_n = eye(n, class(X));
        K   = X*X.' + lambda*I_n;
        % Cholesky solve for K * Alpha = y
        try
            R = chol(K);
            Alpha = R \ (R.' \ y);
        catch
            % SVD fallback (numerically safest)
            [U,S,V] = svd(K, 'econ');
            s = diag(S);
            Alpha = V * ((U.'*y) ./ max(s, eps(class(X))));
        end
        W = X.' * Alpha;
        return;
    end
    % Non-uniform reg requires p×p path
end

% p <= n OR non-uniform reg_mask: use normal equations with Cholesky
A = X.'*X + lambda * diag(reg_mask);
B = X.'*y;

% Symmetrize A slightly to guard against tiny asymmetries
A = (A + A.') * 0.5;

try
    R = chol(A);        % A should be SPD when lambda>0 on all penalized dims
    W = R \ (R.' \ B);
catch
    % SVD fallback: ridge via spectral filter
    % X = U*S*V', ridge solution W = V * diag(s./(s.^2 + lambda_eff)) * U' * y
    [U,S,V] = svd(X, 'econ');
    s = diag(S);
    % If reg_mask not all ones, approximate by mapping to V-basis and scaling:
    % lambda_eff per component ~= lambda * mean(reg_mask)  (simple, conservative)
    if all(abs(reg_mask - 1) < eps)
        coeff = s ./ (s.^2 + lambda);
        W = V * (diag(coeff) * (U.' * y));
    else
        lambda_eff = lambda * mean(reg_mask(:));
        coeff = s ./ (s.^2 + lambda_eff);
        W = V * (diag(coeff) * (U.' * y));
        warning('ridge_regress:svdMaskApprox', ...
            'Used SVD fallback with approximate handling of non-uniform reg_mask.');
    end
end
end
