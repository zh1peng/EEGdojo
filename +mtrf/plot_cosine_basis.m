function plot_cosine_basis(B, lags_ms, varargin)
% PLOT_COSINE_BASIS  Visualizes the set of raised cosine basis functions.
%
%   plot_cosine_basis(B, lags_ms)
%   plot_cosine_basis(..., 'ax', ax_handle)
%   plot_cosine_basis(..., 'title', 'My Title')
%
% Inputs:
%   B         - [L x K] matrix of K basis functions over L lags.
%   lags_ms   - [L x 1] vector of time lags in milliseconds.
%
% Name-Value Pairs:
%   'ax'      - Axes handle to plot into. If empty, a new figure is created.
%   'title'   - Title for the plot.
%
% See also: rcb.build_cosine_basis

ip = inputParser;
ip.addParameter('ax', [], @(h) isempty(h) || isgraphics(h, 'axes'));
ip.addParameter('title', 'Raised-Cosine Basis Functions', @ischar);
ip.parse(varargin{:});
opt = ip.Results;

if isempty(opt.ax)
    figure('Color', 'w');
    ax = axes;
else
    ax = opt.ax;
end

plot(ax, lags_ms(:), B, 'LineWidth', 1.5);
box(ax, 'on');
xlabel(ax, 'Lag (ms)');
ylabel(ax, 'Amplitude');
title(ax, opt.title);
grid(ax, 'on');

end
