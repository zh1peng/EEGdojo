function plot_sim1b
% colored (pink)  Gaussian signal and noise
% varying K
%
% Stefan Haufe, 2017

load('mat/sim1b_results')


figure(1); clf
h0 = plot([min(Ks) max(Ks)], [min(Ks) max(Ks)], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
hold on
set(gca, 'ColorOrderIndex', 1)
h1 = plot(Ks, [mean(K_shift, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h1, 'MarkerFaceColor', get(h1, 'Color'));
patch([Ks'; fliplr(Ks)'], [mean(K_shift, 2) + std(K_shift, [], 2); flipud(mean(K_shift, 2) - std(K_shift, [], 2))], ...
  get(h1, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h2 = plot(Ks, [mean(K_phase, 2) ], 's-', 'linewidth', 2, 'Markersize', 8);
set(h2, 'MarkerFaceColor', get(h2, 'Color'));
patch([Ks'; fliplr(Ks)'], [mean(K_phase, 2) + std(K_phase, [], 2); flipud(mean(K_phase, 2) - std(K_phase, [], 2))], ...
  get(h2, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h3 = plot(Ks, [mean(K_param, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h3, 'MarkerFaceColor', get(h3, 'Color'));
patch([Ks'; fliplr(Ks)'], [mean(K_param, 2) + std(K_param, [], 2); flipud(mean(K_param, 2) - std(K_param, [], 2))], ...
  get(h3, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
set(gca, 'fontsize', 20, 'xtick', [1 5 10 15 20 25])
grid on
axis equal
xlim([min(Ks) max(Ks)])
ylim([min(Ks) max(Ks)])
xlabel('True K')
ylabel('Estimated K')
% title({'Pink Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', N=' num2str(N) ', SNR=' num2str(db_snr) 'dB']})
title('Pink')
legend([h1 h2 h3 h0], 'Shift', 'Phase', 'Parametric', 'Ideal', 'Location', 'Northwest')
ax = reshape(axis, 2, 2);
text(ax(1, 1) - diff(ax(:, 1))*0.125, ax(1, 2) + diff(ax(:, 2))*1.025, 'C', 'fontsize', 30)
export_fig('../../figures/simulations/sim1b_K', '-r300', '-a2');

figure(2); clf
h0 = plot([min(Ks) max(Ks)], [1 1], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
hold on
set(gca, 'ColorOrderIndex', 4)
h1 = plot(Ks, [mean(ISC_av, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h1, 'MarkerFaceColor', get(h1, 'Color'));
patch([Ks'; fliplr(Ks)'], [mean(ISC_av, 2) + std(ISC_av, [], 2); flipud(mean(ISC_av, 2) - std(ISC_av, [], 2))], ...
  get(h1, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h4 = plot(Ks, [mean(ISC_av_test, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h4, 'MarkerFaceColor', get(h4, 'Color'));
patch([Ks'; fliplr(Ks)'], [mean(ISC_av_test, 2) + std(ISC_av_test, [], 2); flipud(mean(ISC_av_test, 2) - std(ISC_av_test, [], 2))], ...
  get(h4, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)

set(gca, 'fontsize', 20)
grid on
xlim([min(Ks) max(Ks)])
ylim([0 1.05])
xlabel('True K')
ylabel('Mean ISC')
% title({'Pink Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', N=' num2str(N) ', SNR=' num2str(db_snr) 'dB']})
title('Pink')
legend([h1 h4 h0], 'Train', 'Test', 'Ideal', 'Location', 'Southwest')
export_fig('../../figures/simulations/sim1b_ISC', '-r300', '-a2');



figure(3); clf
h0 = plot([min(Ks) max(Ks)], [1 1], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
hold on
set(gca, 'ColorOrderIndex', 6)
h2 = plot(Ks, [mean(Y_acc, 2) ], 's-', 'linewidth', 2, 'Markersize', 8);
set(h2, 'MarkerFaceColor', get(h2, 'Color'));
patch([Ks'; fliplr(Ks)'], [mean(Y_acc, 2) + std(Y_acc, [], 2); flipud(mean(Y_acc, 2) - std(Y_acc, [], 2))], ...
  get(h2, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h3 = plot(Ks, [mean(A_acc, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h3, 'MarkerFaceColor', get(h3, 'Color'));
patch([Ks'; fliplr(Ks)'], [mean(A_acc, 2) + std(A_acc, [], 2); flipud(mean(A_acc, 2) - std(A_acc, [], 2))], ...
  get(h3, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
set(gca, 'ColorOrderIndex', 6)
h5 = plot(Ks, [mean(Y_subspace, 2) ], 'd--', 'linewidth', 2, 'Markersize', 8);
set(h5, 'MarkerFaceColor', get(h5, 'Color'));
patch([Ks'; fliplr(Ks)'], [mean(Y_subspace, 2) + std(Y_subspace, [], 2); flipud(mean(Y_subspace, 2) - std(Y_subspace, [], 2))], ...
  get(h5, 'Color'), 'edgecolor', 'none', 'facealpha', 0.05)
h6 = plot(Ks, [mean(A_subspace, 2)], 'd--', 'linewidth', 2, 'Markersize', 8);
set(h6, 'MarkerFaceColor', get(h6, 'Color'));
patch([Ks'; fliplr(Ks)'], [mean(A_subspace, 2) + std(A_subspace, [], 2); flipud(mean(A_subspace, 2) - std(A_subspace, [], 2))], ...
  get(h6, 'Color'), 'edgecolor', 'none', 'facealpha', 0.05)
set(gca, 'fontsize', 20)
grid on
xlim([min(Ks) max(Ks)])
ylim([0 1.05])
xlabel('True K')
ylabel('Reconstruction performance')
% title({'Pink Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', N=' num2str(N) ', SNR=' num2str(db_snr) 'dB']})
title('Pink')
legend([h2 h3 h5 h6 h0], 'Y (component)', 'A (component)', 'Y (subspace)', 'A (subspace)', 'Ideal', 'Location', 'Northeast')
export_fig('../../figures/simulations/sim1b_recon', '-r300', '-a2');

