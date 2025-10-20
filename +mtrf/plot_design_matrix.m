function plot_design_matrix(Xb, F, K, varargin)
% PLOT_DESIGN_MATRIX Visualizes the design matrix using imagesc.
%
%   plot_design_matrix(Xb, F, K)
%   plot_design_matrix(..., 'ax', ax_handle)
%   plot_design_matrix(..., 'title', 'My Title')
%
% Inputs:
%   Xb - [T x (F*K)] design matrix.
%   F  - Number of features.
%   K  - Number of basis functions per feature.
%
% Name-Value Pairs:
%   'ax'    - Axes handle to plot into. If empty, a new figure is created.
%   'title' - Title for the plot.

ip = inputParser;
ip.addParameter('ax', [], @(h) isempty(h) || isgraphics(h, 'axes'));
ip.addParameter('title', 'Design Matrix (Feature x Basis)', @ischar);
ip.parse(varargin{:});
opt = ip.Results;

if isempty(opt.ax)
    figure('Color', 'w');
    ax = axes;
else
    ax = opt.ax;
end

imagesc(ax, Xb);
axis(ax, 'tight');
colorbar(ax);
xlabel(ax, sprintf('Predictors (F=%d features x K=%d bases)', F, K));
ylabel(ax, 'Time (samples)');
title(ax, opt.title);

% Add vertical lines to separate features
hold(ax, 'on');
for f = 1:F-1
    xline(ax, f * K + 0.5, ':k', 'HandleVisibility', 'off');
end
hold(ax, 'off');

end
