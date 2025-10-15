function W_trf = reconstruct_trf(W_basis, B, F, K)
% RECONSTRUCT_TRF Maps basis weights back to the full TRF in lag space.
%
%   W_trf = reconstruct_trf(W_basis, B, F, K)
%
%   The model is fit on weights for the basis functions (W_basis). This
%   function reconstructs the full temporal response function (TRF) for each
%   feature by taking a weighted sum of the basis functions.
%
% Inputs:
%   W_basis - [(F*K) x 1] vector of weights from the regression model.
%   B       - [L x K] basis set matrix.
%   F       - Number of features.
%   K       - Number of basis functions.
%
% Output:
%   W_trf   - [L x F] matrix, where each column is the reconstructed TRF
%             for a single feature.

% Reshape the flat weight vector into a [K x F] matrix,
% where each column contains the K basis weights for one feature.
w_fk = reshape(W_basis, K, F);

% Project the basis weights back into the lag domain.
% B [L x K] * w_fk [K x F] results in W_trf [L x F].
W_trf = B * w_fk;

end