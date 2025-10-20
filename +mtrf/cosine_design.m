function Xb = cosine_design(Xraw, B)
% CONVOLVE_BASIS_DESIGN Convolves features with basis functions.
%
%   Xb = convolve_basis_design(Xraw, B)
%
%   This function creates the design matrix for a TRF model using a basis set.
%   Instead of creating a time-lagged matrix and then projecting it onto the
%   basis, it directly convolves each raw feature channel with each basis
%   function. This is a more memory-efficient approach.
%
% Inputs:
%   Xraw - [T x F] matrix of raw features (T samples, F features).
%   B    - [L x K] basis set, where L is the number of lags and K is the
%          number of basis functions.
%
% Output:
%   Xb   - [T x (F*K)] design matrix, where each feature is represented by
%          its projection onto the K basis functions.

[T, F] = size(Xraw);
K = size(B, 2);

% Pre-allocate the design matrix
Xb = zeros(T, F * K, 'like', Xraw);

% Loop over features and basis functions
for f = 1:F
    x = Xraw(:, f);
    for k = 1:K
        % The convolution of a feature with a basis function is equivalent
        % to the dot product of that basis function with the time-lagged
        % feature matrix at each time point.
        % flipud(B(:,k)) is used to ensure the convolution reflects a
        % causal relationship (i.e., modeling response based on past events).
        
        col_idx = (f - 1) * K + k;
        Xb(:, col_idx) = conv(x, flipud(B(:, k)), 'same');
    end
end

end