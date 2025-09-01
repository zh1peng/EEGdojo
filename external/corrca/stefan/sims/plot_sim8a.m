function plot_sim8a
% Gaussian white signal and noise time series,
% spatial correlation structure of signal/noise either the same or 
% different across subjects
%
% Stefan Haufe, 2017

load('mat/sim8a_results')

xx = (1:4);

cols = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


figure(1); clf
h0 = plot([0.5 4.5], [K K], '--', 'linewidth', 2, 'Color', 'k');
hold on
co = cols(1, :);
for ii = 1:4
  h1 = plot([-0.33; -0.11] + xx(ii), repmat(mean(K_shift(ii, :), 2), 2, 1), '-', 'linewidth', 2);
  set(h1, 'Color', co);
  r = rectangle('Position', [xx(ii)-0.33, mean(K_shift(ii, :), 2)-std(K_shift(ii, :), [], 2), 0.22, 2*std(K_shift(ii, :), [], 2)]);
  set(r, 'FaceColor', [co 0.1], 'EdgeColor', 'none');
end
co = cols(2, :);
for ii = 1:4
  h2 = plot([-0.11; 0.11] + xx(ii), repmat(mean(K_phase(ii, :), 2), 2, 1), '-', 'linewidth', 2);
  set(h2, 'Color', co);
  r = rectangle('Position', [xx(ii)-0.11, mean(K_phase(ii, :), 2)-std(K_phase(ii, :), [], 2), 0.22, 2*std(K_phase(ii, :), [], 2)]);
  set(r, 'FaceColor', [co 0.1], 'EdgeColor', 'none');
end
co = cols(3, :);
for ii = 1:4
  h3 = plot([0.11; 0.33] + xx(ii), repmat(mean(K_param(ii, :), 2), 2, 1), '-', 'linewidth', 2);
  set(h3, 'Color', co);
  r = rectangle('Position', [xx(ii)+0.11, mean(K_param(ii, :), 2)-std(K_param(ii, :), [], 2), 0.22, 2*std(K_param(ii, :), [], 2)]);
  set(r, 'FaceColor', [co 0.1], 'EdgeColor', 'none');
end

set(gca, 'fontsize', 14)
grid on
xlim([0.5 4.5])
ylim([0 13])
set(gca, 'xtick', 1:4, 'xticklabel', {'same/same', 'same/diff.', 'diff./same', 'diff./diff.'})
xlabel('Signal/noise subspaces across subjects', 'FontSize', 20)
ylabel('Estimated K', 'FontSize', 20)
% title({'IID Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', N=' num2str(N) ', K=' num2str(K) ', SNR=' num2str(db_snr) 'dB']}, 'FontSize', 20)
title(['IID, N=' num2str(N)])
l = legend([h1 h2 h3 h0], 'Shift', 'Phase', 'Parametric', 'Ideal', 'Location', 'Southwest');
set(l, 'FontSize', 20)
ax = reshape(axis, 2, 2);
text(ax(1, 1) - diff(ax(:, 1))*0.125, ax(1, 2) + diff(ax(:, 2))*1.025, 'C', 'fontsize', 30)
export_fig('../../figures/simulations/sim8a_K', '-r300', '-a2');




figure(2); clf
h0 = plot([0.5 4.5], [1 1], '--', 'linewidth', 2, 'Color', 'k');
hold on
co = cols(4, :);
for ii = 1:4
  h1 = plot([-0.33; 0] + xx(ii), repmat(nanmean(ISC_av(ii, :), 2), 2, 1), '-', 'linewidth', 2);
  set(h1, 'Color', co);
  r = rectangle('Position', [xx(ii)-0.33, nanmean(ISC_av(ii, :), 2)-nanstd(ISC_av(ii, :), [], 2), 0.33, 2*nanstd(ISC_av(ii, :), [], 2)]);
  set(r, 'FaceColor', [co 0.1], 'EdgeColor', 'none');
end
co = cols(5, :);
for ii = 1:4
  h3 = plot([0; 0.33] + xx(ii), repmat(nanmean(ISC_av_test(ii, :), 2), 2, 1), '-', 'linewidth', 2);
  set(h3, 'Color', co);
  r = rectangle('Position', [xx(ii), nanmean(ISC_av_test(ii, :), 2)-nanstd(ISC_av_test(ii, :), [], 2), 0.33, 2*nanstd(ISC_av_test(ii, :), [], 2)]);
  set(r, 'FaceColor', [co 0.1], 'EdgeColor', 'none');
end

set(gca, 'fontsize', 14)
grid on
xlim([0.5 4.5])
ylim([0 1.05])
set(gca, 'xtick', 1:4, 'xticklabel', {'same/same', 'same/diff.', 'diff./same', 'diff./diff.'})
xlabel('Signal/noise subspaces across subjects', 'FontSize', 20)
ylabel('Mean ISC', 'FontSize', 20)
% title({'IID Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', N=' num2str(N) ', K=' num2str(K) ', SNR=' num2str(db_snr) 'dB']}, 'FontSize', 20)
title(['IID, N=' num2str(N)])
l = legend([h1 h3 h0], 'Train', 'Test', 'Ideal', 'Location', 'Southwest');
set(l, 'FontSize', 20)
ax = reshape(axis, 2, 2);
text(ax(1, 1) - diff(ax(:, 1))*0.125, ax(1, 2) + diff(ax(:, 2))*1.025, 'D', 'fontsize', 30)
export_fig('../../figures/simulations/sim8a_ISC', '-r300', '-a2');





figure(3); clf
h0 = plot([0.5 4.5], [1 1], '--', 'linewidth', 2, 'Color', 'k');
hold on
co = cols(6, :);
for ii = 1:4
  h1 = plot([-0.33; 0] + xx(ii), repmat(nanmean(Y_acc(ii, :), 2), 2, 1), '-', 'linewidth', 2);
  set(h1, 'Color', co);
  r = rectangle('Position', [xx(ii)-0.33, nanmean(Y_acc(ii, :), 2)-nanstd(Y_acc(ii, :), [], 2), 0.33, 2*nanstd(Y_acc(ii, :), [], 2)]);
  set(r, 'FaceColor', [co 0.1], 'EdgeColor', 'none');
end
co = cols(6, :);
for ii = 1:4
  h2 = plot([0; 0.33] + xx(ii), repmat(nanmean(Y_subspace(ii, :), 2), 2, 1), ':', 'linewidth', 2);
  set(h2, 'Color', co);
  r = rectangle('Position', [xx(ii), nanmean(Y_subspace(ii, :), 2)-nanstd(Y_subspace(ii, :), [], 2), 0.33, 2*nanstd(Y_subspace(ii, :), [], 2)]);
  set(r, 'FaceColor', [co 0.05], 'EdgeColor', 'none');
end

set(gca, 'fontsize', 14)
grid on
xlim([0.5 4.5])
ylim([0 1.05])
set(gca, 'xtick', 1:4, 'xticklabel', {'same/same', 'same/diff.', 'diff./same', 'diff./diff.'})
xlabel('Signal/noise subspaces across subjects', 'FontSize', 20)
ylabel('Reconstruction performance', 'FontSize', 20)
% title({'IID Gaussian signal and noise,', ...
%   ['T=' num2str(T) ', D=' num2str(D) ', N=' num2str(N) ', K=' num2str(K) ', SNR=' num2str(db_snr) 'dB']}, 'FontSize', 20)
title(['IID, N=' num2str(N)])
l = legend([h1 h2 h0], 'Y (component)', 'Y (subspace)', 'Ideal', 'Location', 'Northeast');
set(l, 'FontSize', 20)
export_fig('../../figures/simulations/sim8a_recon', '-r300', '-a2');




