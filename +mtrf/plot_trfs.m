function plot_trfs(W_trf, lags_ms, featNames, varargin)
% PLOT_TRFS  Visualizes the reconstructed temporal response functions (TRFs).
%
%   plot_trfs(W_trf, lags_ms, featNames)
%   plot_trfs(..., 'ax', ax_handle)
%   plot_trfs(..., 'title', 'My Title')
%
% Inputs:
%   W_trf     - [L x F] matrix of reconstructed TRFs.
%   lags_ms   - [L x 1] vector of time lags in milliseconds.
%   featNames - [1 x F] cell array or string array of feature names.
%
% Name-Value Pairs:
%   'ax'      - Axes handle to plot into. If empty, a new figure is created.
%   'title'   - Title for the plot.
%
% See also: rcb.reconstruct_trf

ip = inputParser;
ip.addParameter('ax', [], @(h) isempty(h) || isgraphics(h, 'axes'));
ip.addParameter('title', 'Reconstructed TRFs', @ischar);
ip.parse(varargin{:});
opt = ip.Results;

if isempty(opt.ax)
    figure('Color', 'w');
    ax = axes;
else
    ax = opt.ax;
end

plot(ax, lags_ms, W_trf, 'LineWidth', 1.5);
box(ax, 'on');
grid(ax, 'on');
hold(ax, 'on');
xline(ax, 0, '--k', 'HandleVisibility', 'off');
yline(ax, 0, '--k', 'HandleVisibility', 'off');
hold(ax, 'off');

xlabel(ax, 'Lag (ms)');
ylabel(ax, 'Amplitude (a.u.)');
title(ax, opt.title);

if ~isempty(featNames)
    legend(ax, featNames, 'Location', 'best');
end

end
