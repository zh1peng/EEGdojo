csvpath='Z:\matlab_toolbox\EEGdojo/example/movie_csv/ThePresent_summary_codes_10Hz_intuitivenames.csv'
T = readtable(csvpath);
assert(any(strcmpi(T.Properties.VariableNames,'seconds')), 'CSV must have a "seconds" column');
t_src = T.seconds(:);

% Candidate feature names (exclude 'seconds')
varNames = T.Properties.VariableNames;

FsTarget=250;
t0=0;
t1=203.074;
 N  = floor((t1 - t0)*FsTarget) + 1;
    t_tgt = (0:N-1)'/FsTarget + t0;
y_src=T.negative(:);
y_tgt = movie.resample_one_feature(t_src, y_src, t_tgt, 'hold');




win = [1 200];          % seconds
ix_src = t_src >= win(1) & t_src <= win(2);
ix_tgt = t_tgt >= win(1) & t_tgt <= win(2);

fprintf('Unique source values in window: '); disp(unique(y_src(ix_src))')
fprintf('Median dt source = %.6f s (should be ~0.1)\n', median(diff(t_src(ix_src))))

% Show any timestamp irregularities and fast flips in the source
[src_dt_hist,edges] = histcounts(diff(t_src(ix_src)),50);
figure('Color','w'); 
subplot(2,1,1); bar(edges(1:end-1), src_dt_hist, 'histc'); grid on;
xlabel('Δt_src (s)'); ylabel('count'); title('Source Δt near 120 s');

subplot(2,1,2); 
stairs(t_src(ix_src), y_src(ix_src), 'LineWidth',1.2); grid on;
xlim(win); xlabel('Time (s)'); ylabel('y\_src'); title('Source series near 120 s');

figure('Color','w');
stairs(t_tgt(ix_tgt), y_tgt(ix_tgt), 'LineWidth',1.2); grid on;
xlim(win); xlabel('Time (s)'); ylabel('y\_tgt'); title('Resampled (250 Hz) near 120 s');






%% ---- Envelope-style regressors aligned to the same window ----
Fs = FsTarget;                  % 250 Hz target rate
dt = 1/Fs;

% Level (0/1/2/...) and binary mask
y_level = y_tgt(:);
y_bin   = y_level > 0;

% 1) Onset stick (rising edges only)
onset = [false; diff(y_bin) > 0];
onset = double(onset);

% 2) Causal exponential envelope: env[t] = alpha*env[t-1] + onset[t]
tau_sec = 0.75;                 % try 0.3–1.5 s
alpha   = exp(-dt/tau_sec);
env_exp = filter(1, [1 -alpha], onset);

% 3) Causal moving-rate (event count in a W-second trailing window)
win_sec = 1.0;                  % trailing window length
W = max(1, round(win_sec * Fs));
kern = ones(W,1);               % causal boxcar
rate_mov = conv(onset, kern, 'same');     % ~events/s * W (not normalized)

% (Optional) normalize envelopes for visualization (within plotted window)
norm_in_win = @(x) x ./ max(x(ix_tgt) + eps);
env_exp_n  = norm_in_win(env_exp);
rate_mov_n = norm_in_win(rate_mov);

%% Figure — Source, resampled level, and envelopes (aligned to `win`)
figure('Color','w','Name','Feature envelopes (aligned)');
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');

% Row 1: source series (10 Hz)
nexttile; hold on;
stairs(t_src(ix_src), y_src(ix_src), 'LineWidth',1.1, 'Color',[0.35 0.35 0.35]);
grid on; xlim(win); ylabel('y_{src}');
title(sprintf('Source series (%.1f–%.1f s)', win(1), win(2)));

% Row 2: resampled levels @ 250 Hz
nexttile; hold on;
stairs(t_tgt(ix_tgt), y_level(ix_tgt), 'LineWidth',1.2, 'Color',[0 0.45 0.74]);
grid on; xlim(win); ylabel('level');
title('Resampled to EEG grid (250 Hz)');

% Row 3: onset stick (rising edges)
nexttile; hold on;
stem(t_tgt(ix_tgt), onset(ix_tgt), 'Marker','none', 'Color',[0.85 0.33 0.10], 'LineWidth',1.0);
grid on; xlim(win); ylim([-0.1 1.1]); ylabel('onset');
title('Onset stick (rising edges)');

% Row 4: causal envelopes (exp + moving rate)
nexttile; hold on;
plot(t_tgt(ix_tgt), env_exp_n(ix_tgt),  'LineWidth',1.4, 'DisplayName',sprintf('exp \\tau=%.2fs', tau_sec));
plot(t_tgt(ix_tgt), rate_mov_n(ix_tgt), '--', 'LineWidth',1.2, 'DisplayName',sprintf('rate %gs win', win_sec));
grid on; xlim(win); ylim([-0.05 1.05]);
xlabel('Time (s)'); ylabel('normalized');
title('Causal “envelopes” (normalized in window)');
legend('Location','best'); legend boxoff;



%% Smooth "contour" versions of the resampled level (y_tgt)
Fs = FsTarget;
lvl = double(y_tgt(:));

% A) causal EMA
tau_sec = 0.6;                         % try 0.3–1.0 s
alpha   = exp(-1/(Fs*tau_sec));
lvl_ema = filter(1, [1 -alpha], lvl);

% B) symmetric Gaussian smoother (viz only)
win_sec   = 0.8;                       % ~window length
win_gauss = max(3, round(win_sec*Fs));
lvl_gauss = smoothdata(lvl, 'gaussian', win_gauss);

% C) Savitzky–Golay (viz only)
sg_win = 2*floor(0.3*Fs)+1;            % ~0.6 s; must be odd
lvl_sg = smoothdata(lvl, 'sgolay', sg_win);

%% Plot in separate subplots, aligned to your window (ix_src, ix_tgt)
figure('Color','w','Name','Smoothed contours vs original/resampled');
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');

% 1) Source vs Resampled
nexttile; hold on;
stairs(t_src(ix_src), y_src(ix_src), 'LineWidth',1.0, 'Color',[0.35 0.35 0.35], ...
    'DisplayName','source (10 Hz)');
stairs(t_tgt(ix_tgt), lvl(ix_tgt),   'LineWidth',1.2, 'Color',[0 0.45 0.74], ...
    'DisplayName','resampled (250 Hz)');
xlim(win); grid on; ylabel('level');
title('Source vs Resampled'); legend('Location','best'); legend boxoff;

% 2) EMA contour (causal) + resampled background
nexttile; hold on;
plot( t_tgt(ix_tgt), lvl(ix_tgt),    'Color',[0.85 0.85 0.85], 'LineWidth',1.0, ...
    'DisplayName','resampled');
plot( t_tgt(ix_tgt), lvl_ema(ix_tgt)/100,'LineWidth',1.6, ...
    'DisplayName',sprintf('EMA  \\tau=%.1fs', tau_sec));
xlim(win); grid on; ylabel('level');
title('Causal EMA contour'); legend('Location','best'); legend boxoff;

% 3) Gaussian contour (non-causal) + resampled background
nexttile; hold on;
plot( t_tgt(ix_tgt), lvl(ix_tgt),      'Color',[0.85 0.85 0.85], 'LineWidth',1.0, ...
    'DisplayName','resampled');
plot( t_tgt(ix_tgt), lvl_gauss(ix_tgt),'--', 'LineWidth',1.6, ...
    'DisplayName',sprintf('Gaussian ~%gs', win_sec));
xlim(win); grid on; ylabel('level');
title('Gaussian-smoothed contour'); legend('Location','best'); legend boxoff;

% 4) Savitzky–Golay contour (non-causal) + resampled background
nexttile; hold on;
plot( t_tgt(ix_tgt), lvl(ix_tgt),   'Color',[0.85 0.85 0.85], 'LineWidth',1.0, ...
    'DisplayName','resampled');
plot( t_tgt(ix_tgt), lvl_sg(ix_tgt), ':', 'LineWidth',1.8, ...
    'DisplayName',sprintf('Savitzky–Golay ~%.1fs', sg_win/Fs));
xlim(win); grid on; ylabel('level'); xlabel('Time (s)');
title('Savitzky–Golay contour'); legend('Location','best'); legend boxoff;
