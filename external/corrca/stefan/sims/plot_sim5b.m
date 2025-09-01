function plot_sim5b
% colored (pink) Gaussian signal and noise
% vary T
%
% Stefan Haufe, 2017

load('mat/sim5b_results')


figure(1); clf
h0 = plot([min(Ns) max(Ns)], [K K], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
hold on
set(gca, 'ColorOrderIndex', 1)
h1 = plot(Ns, [mean(K_shift, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h1, 'MarkerFaceColor', get(h1, 'Color'));
patch([Ns'; fliplr(Ns)'], [mean(K_shift, 2) + std(K_shift, [], 2); flipud(mean(K_shift, 2) - std(K_shift, [], 2))], ...
  get(h1, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h2 = plot(Ns, [mean(K_phase, 2) ], 's-', 'linewidth', 2, 'Markersize', 8);
set(h2, 'MarkerFaceColor', get(h2, 'Color'));
patch([Ns'; fliplr(Ns)'], [mean(K_phase, 2) + std(K_phase, [], 2); flipud(mean(K_phase, 2) - std(K_phase, [], 2))], ...
  get(h2, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h3 = plot(Ns, [mean(K_param, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h3, 'MarkerFaceColor', get(h3, 'Color'));
patch([Ns'; fliplr(Ns)'], [mean(K_param, 2) + std(K_param, [], 2); flipud(mean(K_param, 2) - std(K_param, [], 2))], ...
  get(h3, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
set(gca, 'Fontsize', 20, 'xtick', [2 5 10 15 20])
grid on
xlim([min(Ns) max(Ns)])
ylim([0 16])
xlabel('N')
ylabel('Estimated K')
% title({'Pink Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', K=' num2str(K) ', SNR=' num2str(db_snr) 'dB']})
title(['Pink, T=' num2str(T) ])
legend([h1 h2 h3 h0], 'Shift', 'Phase', 'Parametric', 'Ideal', 'Location', 'Southeast')
ax = reshape(axis, 2, 2);
text(ax(1, 1) - diff(ax(:, 1))*0.125, ax(1, 2) + diff(ax(:, 2))*1.025, 'D', 'fontsize', 30)
export_fig('../../figures/simulations/sim5b_K', '-r300', '-a2');

figure(2); clf
h0 = plot([min(Ns) max(Ns)], [1 1], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
hold on
set(gca, 'ColorOrderIndex', 4)
h1 = plot(Ns, [mean(ISC_av, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h1, 'MarkerFaceColor', get(h1, 'Color'));
patch([Ns'; fliplr(Ns)'], [mean(ISC_av, 2) + std(ISC_av, [], 2); flipud(mean(ISC_av, 2) - std(ISC_av, [], 2))], ...
  get(h1, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h4 = plot(Ns, [mean(ISC_av_test, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h4, 'MarkerFaceColor', get(h4, 'Color'));
patch([Ns'; fliplr(Ns)'], [mean(ISC_av_test, 2) + std(ISC_av_test, [], 2); flipud(mean(ISC_av_test, 2) - std(ISC_av_test, [], 2))], ...
  get(h4, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)

set(gca, 'Fontsize', 20, 'xtick', [2 5 10 15 20])
grid on
xlim([min(Ns) max(Ns)])
ylim([0 1.05])
xlabel('N')
ylabel('Mean ISC')
% title({'Pink Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', K=' num2str(K) ', SNR=' num2str(db_snr) 'dB']})
title(['Pink, T=' num2str(T) ])
legend([h1 h4 h0], 'Train', 'Test', 'Ideal', 'Location', 'Southeast')
export_fig('../../figures/simulations/sim5b_ISC', '-r300', '-a2');



figure(3); clf
h0 = plot([min(Ns) max(Ns)], [1 1], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
hold on
set(gca, 'ColorOrderIndex', 6)
h2 = plot(Ns, [mean(Y_acc, 2) ], 's-', 'linewidth', 2, 'Markersize', 8);
set(h2, 'MarkerFaceColor', get(h2, 'Color'));
patch([Ns'; fliplr(Ns)'], [mean(Y_acc, 2) + std(Y_acc, [], 2); flipud(mean(Y_acc, 2) - std(Y_acc, [], 2))], ...
  get(h2, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
h3 = plot(Ns, [mean(A_acc, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
set(h3, 'MarkerFaceColor', get(h3, 'Color'));
patch([Ns'; fliplr(Ns)'], [mean(A_acc, 2) + std(A_acc, [], 2); flipud(mean(A_acc, 2) - std(A_acc, [], 2))], ...
  get(h3, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
set(gca, 'ColorOrderIndex', 6)
h5 = plot(Ns, [mean(Y_subspace, 2) ], 'd--', 'linewidth', 2, 'Markersize', 8);
set(h5, 'MarkerFaceColor', get(h5, 'Color'));
patch([Ns'; fliplr(Ns)'], [mean(Y_subspace, 2) + std(Y_subspace, [], 2); flipud(mean(Y_subspace, 2) - std(Y_subspace, [], 2))], ...
  get(h5, 'Color'), 'edgecolor', 'none', 'facealpha', 0.05)
h6 = plot(Ns, [mean(A_subspace, 2)], 'd--', 'linewidth', 2, 'Markersize', 8);
set(h6, 'MarkerFaceColor', get(h6, 'Color'));
patch([Ns'; fliplr(Ns)'], [mean(A_subspace, 2) + std(A_subspace, [], 2); flipud(mean(A_subspace, 2) - std(A_subspace, [], 2))], ...
  get(h6, 'Color'), 'edgecolor', 'none', 'facealpha', 0.05)
set(gca, 'Fontsize', 20, 'xtick', [2 5 10 15 20])
grid on
xlim([min(Ns) max(Ns)])
ylim([0 1.05])
xlabel('N')
ylabel('Reconstruction performance')
% title({'Pink Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', K=' num2str(K) ', SNR=' num2str(db_snr) 'dB']})
title(['Pink, T=' num2str(T) ])
legend([h2 h3 h5 h6 h0], 'Y (component)', 'A (component)', 'Y (subspace)', 'A (subspace)', 'Ideal', 'Location', 'Southeast')
export_fig('../../figures/simulations/sim5b_recon', '-r300', '-a2');

