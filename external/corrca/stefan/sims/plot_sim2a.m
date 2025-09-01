function plot_sim2a
% white Gaussian signal and noise
% varying snr
%
% Stefan Haufe, 2017

load('mat/sim2a_results')

% co = [         0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840];

figure(1); clf
h0 = plot([min(db_snrs) max(db_snrs)], [K K], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
hold on
set(gca, 'ColorOrderIndex', 1)
h1 = plot(db_snrs, [mean(K_shift, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h1, 'MarkerFaceColor', get(h1, 'Color'));
patch([db_snrs'; fliplr(db_snrs)'], [mean(K_shift, 2) + std(K_shift, [], 2); flipud(mean(K_shift, 2) - std(K_shift, [], 2))], ...
  get(h1, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h2 = plot(db_snrs, [mean(K_phase, 2) ], 's-', 'linewidth', 2, 'Markersize', 8);
set(h2, 'MarkerFaceColor', get(h2, 'Color'));
patch([db_snrs'; fliplr(db_snrs)'], [mean(K_phase, 2) + std(K_phase, [], 2); flipud(mean(K_phase, 2) - std(K_phase, [], 2))], ...
  get(h2, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h3 = plot(db_snrs, [mean(K_param, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h3, 'MarkerFaceColor', get(h3, 'Color'));
patch([db_snrs'; fliplr(db_snrs)'], [mean(K_param, 2) + std(K_param, [], 2); flipud(mean(K_param, 2) - std(K_param, [], 2))], ...
  get(h3, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
set(gca, 'fontsize', 20)
grid on
xlim([min(db_snrs) max(db_snrs)])
ylim([0 13])
xlabel('SNR [dB]')
ylabel('Estimated K')
% title({'IID Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', N=' num2str(N) ', K=' num2str(K)]})
title('IID')
legend([h1 h2 h3 h0], 'Shift', 'Phase', 'Parametric', 'Ideal', 'Location', 'Southeast')
ax = reshape(axis, 2, 2);
text(ax(1, 1) - diff(ax(:, 1))*0.125, ax(1, 2) + diff(ax(:, 2))*1.025, 'C', 'fontsize', 30)
export_fig('../../figures/simulations/sim2a_K', '-r300', '-a2');

print('../../figures/simulations/sim2a_K_test', '-depsc');

figure(2); clf
h0 = plot([min(db_snrs) max(db_snrs)], [1 1], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
hold on
set(gca, 'ColorOrderIndex', 4)
h1 = plot(db_snrs, [mean(ISC_av, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h1, 'MarkerFaceColor', get(h1, 'Color'));
patch([db_snrs'; fliplr(db_snrs)'], [mean(ISC_av, 2) + std(ISC_av, [], 2); flipud(mean(ISC_av, 2) - std(ISC_av, [], 2))], ...
  get(h1, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h4 = plot(db_snrs, [mean(ISC_av_test, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h4, 'MarkerFaceColor', get(h4, 'Color'));
patch([db_snrs'; fliplr(db_snrs)'], [mean(ISC_av_test, 2) + std(ISC_av_test, [], 2); flipud(mean(ISC_av_test, 2) - std(ISC_av_test, [], 2))], ...
  get(h4, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)

set(gca, 'fontsize', 20)
grid on
xlim([min(db_snrs) max(db_snrs)])
ylim([0 1.05])
xlabel('SNR [dB]')
ylabel('Mean ISC')
% title({'IID Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', N=' num2str(N) ', K=' num2str(K)]})
title('IID')
legend([h1 h4 h0], 'Train', 'Test', 'Ideal', 'Location', 'Southeast')
ax = reshape(axis, 2, 2);
text(ax(1, 1) - diff(ax(:, 1))*0.125, ax(1, 2) + diff(ax(:, 2))*1.025, 'A', 'fontsize', 30)
export_fig('../../figures/simulations/sim2a_ISC', '-r300', '-a2');



figure(3); clf
h0 = plot([min(db_snrs) max(db_snrs)], [1 1], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
hold on
set(gca, 'ColorOrderIndex', 6)
h2 = plot(db_snrs, [mean(Y_acc, 2) ], 's-', 'linewidth', 2, 'Markersize', 8);
set(h2, 'MarkerFaceColor', get(h2, 'Color'));
patch([db_snrs'; fliplr(db_snrs)'], [mean(Y_acc, 2) + std(Y_acc, [], 2); flipud(mean(Y_acc, 2) - std(Y_acc, [], 2))], ...
  get(h2, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h3 = plot(db_snrs, [mean(A_acc, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h3, 'MarkerFaceColor', get(h3, 'Color'));
patch([db_snrs'; fliplr(db_snrs)'], [mean(A_acc, 2) + std(A_acc, [], 2); flipud(mean(A_acc, 2) - std(A_acc, [], 2))], ...
  get(h3, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
set(gca, 'ColorOrderIndex', 6)
h5 = plot(db_snrs, [mean(Y_subspace, 2) ], 'd--', 'linewidth', 2, 'Markersize', 8);
set(h5, 'MarkerFaceColor', get(h5, 'Color'));
patch([db_snrs'; fliplr(db_snrs)'], [mean(Y_subspace, 2) + std(Y_subspace, [], 2); flipud(mean(Y_subspace, 2) - std(Y_subspace, [], 2))], ...
  get(h5, 'Color'), 'edgecolor', 'none', 'facealpha', 0.05)
h6 = plot(db_snrs, [mean(A_subspace, 2)], 'd--', 'linewidth', 2, 'Markersize', 8);
set(h6, 'MarkerFaceColor', get(h6, 'Color'));
patch([db_snrs'; fliplr(db_snrs)'], [mean(A_subspace, 2) + std(A_subspace, [], 2); flipud(mean(A_subspace, 2) - std(A_subspace, [], 2))], ...
  get(h6, 'Color'), 'edgecolor', 'none', 'facealpha', 0.05)
set(gca, 'fontsize', 20)
grid on
xlim([min(db_snrs) max(db_snrs)])
ylim([0 1.05])
xlabel('SNR [dB]')
ylabel('Reconstruction performance')
% title({'IID Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', N=' num2str(N) ', K=' num2str(K)]})
title('IID')
legend([h2 h3 h5 h6 h0], 'Y (component)', 'A (component)', 'Y (subspace)', 'A (subspace)', 'Ideal', 'Location', 'Southeast')
ax = reshape(axis, 2, 2);
text(ax(1, 1) - diff(ax(:, 1))*0.125, ax(1, 2) + diff(ax(:, 2))*1.025, 'B', 'fontsize', 30)
export_fig('../../figures/simulations/sim2a_recon', '-r300', '-a2');
