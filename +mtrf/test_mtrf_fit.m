%% ===== mtrf_fit quick end-to-end tests (paste into Command Window) =====
rng(42);

% --- dispatch: use mtrf.mtrf_fit if available, else mtrf_fit

mfit = @(Y,X,fs,varargin) mtrf.mtrf_fit(Y,X,fs,varargin{:});

   


% -------------------------- design ----------------------------
fs       = 100;                     % Hz
lags_ms  = -100:10:500;             % TRF window (ms)
L        = numel(lags_ms);
T        = 6000;                    % 60 s at 100 Hz
F        = 2;                       % features
M        = 1;                       % neural channels

% --------------------- ground-truth TRF -----------------------
biph = @(L,mu,w,a1,a2) (a1*exp(-0.5*((L-mu)/max(w,1)).^2) + a2*exp(-0.5*((L-(mu+80))/max(w*1.4,1)).^2));
h1 = biph(lags_ms(:), +120,  60, +1.0, -0.6);  h1 = h1 ./ max(abs(h1)+eps);
h2 = biph(lags_ms(:), +180, 100, +0.7, -0.4);  h2 = 0.6 * (h2 ./ max(abs(h2)+eps));
Htrue = [h1, h2];                   % [L x F]

% ---------------- stimulus & response -------------------------
Xraw = 0.5*randn(T,F);
Xraw = filter(1, [1 -0.7], Xraw);   % AR(1)-ish color

% Lagged convolution (explicit, no functions) to create noise-free Y
lags_s = round(lags_ms * fs / 1000);
Ysig = zeros(T,M);
for ii = 1:numel(lags_s)
    lag = lags_s(ii);
    if lag>=0
        if lag < T
            Ysig(1+lag:end,1) = Ysig(1+lag:end,1) ...
                + Htrue(ii,1) * Xraw(1:end-lag,1) ...
                + Htrue(ii,2) * Xraw(1:end-lag,2);
        end
    else
        lag2 = -lag;
        if lag2 < T
            Ysig(1:end-lag2,1) = Ysig(1:end-lag2,1) ...
                + Htrue(ii,1) * Xraw(1+lag2:end,1) ...
                + Htrue(ii,2) * Xraw(1+lag2:end,2);
        end
    end
end

% Two SNR settings (per-output scaling by std of Ysig)
snr_good = 6;   % easier
snr_bad  = 0.8; % harder
noise_g  = randn(T,M) * (std(Ysig,0,1) / max(snr_good,eps));
noise_b  = randn(T,M) * (std(Ysig,0,1) / max(snr_bad, eps));
Y_good   = Ysig + noise_g;
Y_bad    = Ysig + noise_b;

% Common CV options
KOuter = 5; KInner = 4;
lambdaGrid = struct('lambda', 2.^(0:2:18));
guardBandSec = max(abs(lags_ms))/1000;   % avoid train/test bleed via lags

% ============================ RUNS ==============================
% 1) Lagged, good SNR
OUT_lag_good = mfit(Y_good, Xraw, fs, ...
    'DesignType','lagged','lags_ms',lags_ms, ...
    'Mode','both','Model','ridge', ...
    'ParamGrid',lambdaGrid, ...
    'SelectMetric','corr', ...
    'KOuter',KOuter,'KInner',KInner, ...
    'GuardBandSec',guardBandSec, ...
    'UseParallel',false,'Verbose',1, ...
    'VariableNames',["feat1","feat2"]);

% 2) Lagged, bad SNR
OUT_lag_bad = mfit(Y_bad, Xraw, fs, ...
    'DesignType','lagged','lags_ms',lags_ms, ...
    'Mode','both','Model','ridge', ...
    'ParamGrid',lambdaGrid, ...
    'SelectMetric','corr', ...
    'KOuter',KOuter,'KInner',KInner, ...
    'GuardBandSec',guardBandSec, ...
    'UseParallel',false,'Verbose',0, ...
    'VariableNames',["feat1","feat2"]);

% 3) Cosine basis on good SNR (smoother TRF)
centers_ms = -100:80:500;  width_ms = 80;
OUT_cosine = mfit(Y_good, Xraw, fs, ...
    'DesignType','cosine', ...
    'lags_ms',lags_ms,'centers_ms',centers_ms,'width_ms',width_ms, ...
    'Mode','both','Model','ridge', ...
    'ParamGrid',lambdaGrid, ...
    'SelectMetric','corr', ...
    'KOuter',KOuter,'KInner',KInner, ...
    'GuardBandSec',guardBandSec, ...
    'UseParallel',false,'Verbose',0, ...
    'VariableNames',["feat1","feat2"]);

% ============================ PLOTS =============================
% ---- TRF recovery (true vs estimated) ----
figure('Name','TRF recovery: lagged vs cosine','Color','w');

subplot(2,2,1); % Lagged good
hold on;
Hhat = OUT_lag_good.deploy.TRF_lag_eff(:,:,1); % [L x F]
plot(lags_ms, Htrue(:,1), '-', 'LineWidth',1.6);
plot(lags_ms, Htrue(:,2), '-', 'LineWidth',1.6);
plot(lags_ms, Hhat(:,1),  '--', 'LineWidth',1.6);
plot(lags_ms, Hhat(:,2),  '--', 'LineWidth',1.6);
grid on; xlabel('Lag (ms)'); ylabel('TRF (orig units)');
title('Lagged (good SNR)'); legend('True f1','True f2','Est f1','Est f2','Location','best');

subplot(2,2,2); % Lagged bad
hold on;
Hhat = OUT_lag_bad.deploy.TRF_lag_eff(:,:,1);
plot(lags_ms, Htrue(:,1), '-', 'LineWidth',1.6);
plot(lags_ms, Htrue(:,2), '-', 'LineWidth',1.6);
plot(lags_ms, Hhat(:,1),  '--', 'LineWidth',1.6);
plot(lags_ms, Hhat(:,2),  '--', 'LineWidth',1.6);
grid on; xlabel('Lag (ms)'); ylabel('TRF (orig units)');
title('Lagged (bad SNR)'); legend('True f1','True f2','Est f1','Est f2','Location','best');

subplot(2,2,3); % Cosine good
hold on;
Hhat = OUT_cosine.deploy.TRF_lag_eff(:,:,1);
plot(lags_ms, Htrue(:,1), '-', 'LineWidth',1.6);
plot(lags_ms, Htrue(:,2), '-', 'LineWidth',1.6);
plot(lags_ms, Hhat(:,1),  '--', 'LineWidth',1.6);
plot(lags_ms, Hhat(:,2),  '--', 'LineWidth',1.6);
grid on; xlabel('Lag (ms)'); ylabel('TRF (orig units)');
title('Cosine basis (good SNR)'); legend('True f1','True f2','Est f1','Est f2','Location','best');

% ---- CV metrics summary (per fold r and text) ----
% ---- CV metrics summary (per fold r and text) ----
figure('Name','CV metrics (mtrf_fit)','Color','w');

% --- Panel 1: Lagged, good SNR ---
subplot(1,3,1);
K = numel(OUT_lag_good.cv.panels.by_fold);
zvals = zeros(1,K);
for kk=1:K
    zvals(kk) = OUT_lag_good.cv.panels.by_fold{kk}.r_fisher_mean;
end
bar(1:K, tanh(zvals)); ylim([-1 1]); grid on;
title('Lagged, good SNR - r by fold'); xlabel('Fold'); ylabel('r');

% --- Panel 2: Lagged, bad SNR ---
subplot(1,3,2);
K = numel(OUT_lag_bad.cv.panels.by_fold);
zvals = zeros(1,K);
for kk=1:K
    zvals(kk) = OUT_lag_bad.cv.panels.by_fold{kk}.r_fisher_mean;
end
bar(1:K, tanh(zvals)); ylim([-1 1]); grid on;
title('Lagged, bad SNR - r by fold'); xlabel('Fold'); ylabel('r');

% --- Panel 3: Cosine, good SNR ---
subplot(1,3,3);
K = numel(OUT_cosine.cv.panels.by_fold);
zvals = zeros(1,K);
for kk=1:K
    zvals(kk) = OUT_cosine.cv.panels.by_fold{kk}.r_fisher_mean;
end
bar(1:K, tanh(zvals)); ylim([-1 1]); grid on;
title('Cosine, good SNR - r by fold'); xlabel('Fold'); ylabel('r');

figure('Name','CV text summaries','Color','w');
subplot(1,3,1); axis off; S = OUT_lag_good.cv.panels.summary;
text(0.05,0.9, {
  sprintf('[Lagged-good]  r:  %.3f ± %.3f', S.r_overall_mean, S.r_overall_se)
  sprintf('                R^2: %.3f ± %.3f', S.R2_mean_mean,  S.R2_mean_se)
  sprintf('                EV:  %.3f ± %.3f', S.EV_mean_mean,  S.EV_mean_se)
}, 'FontName','Consolas','VerticalAlignment','top');

subplot(1,3,2); axis off; S = OUT_lag_bad.cv.panels.summary;
text(0.05,0.9, {
  sprintf('[Lagged-bad]   r:  %.3f ± %.3f', S.r_overall_mean, S.r_overall_se)
  sprintf('                R^2: %.3f ± %.3f', S.R2_mean_mean,  S.R2_mean_se)
  sprintf('                EV:  %.3f ± %.3f', S.EV_mean_mean,  S.EV_mean_se)
}, 'FontName','Consolas','VerticalAlignment','top');

subplot(1,3,3); axis off; S = OUT_cosine.cv.panels.summary;
text(0.05,0.9, {
  sprintf('[Cosine-good]  r:  %.3f ± %.3f', S.r_overall_mean, S.r_overall_se)
  sprintf('                R^2: %.3f ± %.3f', S.R2_mean_mean,  S.R2_mean_se)
  sprintf('                EV:  %.3f ± %.3f', S.EV_mean_mean,  S.EV_mean_se)
}, 'FontName','Consolas','VerticalAlignment','top');

% ---- OUTER-TEST overlay (z-space) ----
figure('Name','OUTER-TEST overlay (lagged, good SNR)','Color','w');
Yz  = OUT_lag_good.cv.y_outer;
Yhz = OUT_lag_good.cv.predict_outer;
mask = ~any(isnan(Yz),2) & ~any(isnan(Yhz),2);
idx  = find(mask);
plot(idx, Yz(mask,1), '-', 'LineWidth', 1); hold on;
plot(idx, Yhz(mask,1), '--', 'LineWidth', 1);
grid on; xlabel('Time (samples)'); ylabel('Z');
legend('Y (z)','Yhat (z)','Location','best');
title('Lagged design: OUTER-TEST, first output');

% ----------------- Console summaries ---------------------------
fprintf('\n[Lagged-good]  '); disp(OUT_lag_good.cv.panels.summary);
fprintf('[Lagged-bad]   '); disp(OUT_lag_bad.cv.panels.summary);
fprintf('[Cosine-good]  '); disp(OUT_cosine.cv.panels.summary);

%% ======================================================================




% Build and orthonormalize basis
B = mtrf.build_cosine(lags_ms, -100:30:500, 60);   % same as above
[Q,~] = qr(B,0);                                   % Q is orthonormal (L x K)

% Cosine design with Q instead of B
Xcos = mtrf.cosine_design(Xraw, Q);                % [T x (F*K)]

% Run nested CV on this custom design (same options as mtrf_fit uses)
OUTQ = mtrf.nested_cv_fit(Y_good, Xcos, ...
  'Mode','both','Model','ridge', ...
  'ParamGrid',struct('lambda',2.^(0:2:18)), ...
  'SelectMetric','corr','KOuter',5,'KInner',4, ...
  'GuardBand',round(max(abs(lags_ms))/1000*fs), ...
  'UseParallel',false,'Verbose',0);

% Map weights back to TRF: for each output m, WfK -> TRF = Q * WfK
Wz   = OUTQ.deploy.model.W;                         % z-space
muX  = OUTQ.deploy.preproc.muX; sdX = OUTQ.deploy.preproc.sdX;
sdY  = OUTQ.deploy.preproc.sdY;                     % (orig-units mapping)
Weff = (Wz .* (1./sdX(:))) .* sdY(:)';              % [(F*K) x M]
Kc   = size(Q,2); F2 = size(Xraw,2); M2 = size(Weff,2);
TRFQ = zeros(numel(lags_ms), F2, M2);
for m=1:M2
    WfK = reshape(Weff(:,m), Kc, F2);               % [K x F]
    TRFQ(:,:,m) = Q * WfK;                          % [L x F]
end

% Plot vs ground truth
figure('Color','w');
plot(lags_ms,Htrue(:,1),'-',lags_ms,Htrue(:,2),'-', ...
     lags_ms,TRFQ(:,1,1),'--',lags_ms,TRFQ(:,2,1),'--','LineWidth',1.6);
grid on; xlabel('Lag (ms)'); ylabel('TRF');
title('Cosine (orthonormal basis)'); legend('True f1','True f2','Est f1','Est f2');
