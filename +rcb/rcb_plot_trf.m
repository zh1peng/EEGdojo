function rcb_plot_trf(out, varargin)
% RCB_PLOT_TRF  Plot TRF waveforms from rcb_fit() output.
% Usage:
%   rcb_plot_trf(out)                       % z-space (default)
%   rcb_plot_trf(out,'Space','original')    % original units
%   rcb_plot_trf(out,'Features',1:6)        % subset of features
%   rcb_plot_trf(out,'Style','image')       % heatmap instead of lines

ip = inputParser;
ip.addParameter('Space','z',@(s) any(strcmpi(s,{'z','original'})));
ip.addParameter('Features',[],@isnumeric);
ip.addParameter('Style','lines',@(s) any(strcmpi(s,{'lines','image'})));
ip.parse(varargin{:});
opt = ip.Results;

B        = out.basis.B;                        % [L x K]
lags_ms  = out.basis.lags_ms;
W_basis  = out.W_basis;                         % [(F*K [+1]) x 1]
featNames= string(out.featNames);

% Handle intercept (if present)
hasIntercept = endsWith(featNames, "intercept");
if any(hasIntercept)
    featNames = featNames(~hasIntercept);
    W_basis   = W_basis(1:end-1,:);
end

% Optional: convert to original units: W_orig = (sd_y / sd_x_j) * W_z
if strcmpi(opt.Space,'original')
    sdY = out.model.sd_y;
    sdX = out.model.sd_X; % vector, length matches columns of Xb used in final fit
    W_basis = (sdY) .* (W_basis ./ sdX(:));
end

% Reconstruct TRF in lag space [L x F]
F = numel(featNames);
K = size(B,2);
W_trf = reshape(W_basis, K, F);        % [K x F]
W_trf = B * W_trf;                     % [L x F]

% Subset features if requested
if ~isempty(opt.Features)
    W_trf     = W_trf(:, opt.Features);
    featNames = featNames(opt.Features);
end

% Plot
switch lower(opt.Style)
    case 'lines'
        figure; hold on;
        for f = 1:numel(featNames)
            plot(lags_ms, W_trf(:,f), 'LineWidth',1.5);
        end
        xline(0,'k:'); yline(0,'k-');
        xlabel('Lag (ms)'); ylabel('TRF amplitude');
        title(sprintf('TRF (%s-space)', opt.Space));
        legend(featNames, 'Interpreter','none','Location','bestoutside');
        grid on;
    case 'image'
        figure;
        imagesc(lags_ms, 1:numel(featNames), W_trf.' );
        axis xy; colorbar;
        xlabel('Lag (ms)'); ylabel('Feature');
        yticks(1:numel(featNames)); yticklabels(featNames);
        title(sprintf('TRF heatmap (%s-space)', opt.Space));
        hold on; xline(0,'k:');
end
end
