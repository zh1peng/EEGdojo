rng(7);
    ncv = @(Y,X,varargin) mtrf.nested_cv_fit(Y,X,varargin{:});

N = 4000;  P = 20;  M = 2;

% Correlated features
SigmaX = toeplitz(0.9.^(0:P-1));
X = randn(N,P) * chol(SigmaX);

% True weights
W = randn(P,M); W = W./vecnorm(W);

Yclean = X*W;

% --- Per-OUTPUT SNR control (better than Frobenius-wide scaling) ---
snr = 6;  % try 6–10 for an easy, obvious effect
noise = randn(N,M);
for j = 1:M
    noise(:,j) = noise(:,j) .* (std(Yclean(:,j)) / snr);
end
Y = Yclean + noise;

% --- Tighter λ-grid (avoid over-shrink) ---
lambdaGrid = struct('lambda', logspace(-2,2,12));  % 0.01 ... 100

OUT = ncv(Y,X,'Mode','assess','Model','ridge', ...
    'ParamGrid', lambdaGrid, ...
    'SelectMetric','corr', ...          % you can also try 'corr'
    'KOuter',5,'KInner',4, ...
    'GuardBand',0,'UseParallel',false,'Verbose',1);

% ---- Quick plots/summary ----
K = numel(OUT.panels.by_fold);
z = zeros(K,1);
for k=1:K, z(k) = OUT.panels.by_fold{k}.r_fisher_mean; end
figure('Color','w'); 
subplot(1,2,1); bar(tanh(z)); ylim([-1 1]); grid on;
xlabel('Outer fold'); ylabel('Overall r'); title('Fisher-mean r (tanh back)');

S = OUT.panels.summary;
subplot(1,2,2); axis off;
text(0.05,0.9, {
  sprintf('r:   %.3f ± %.3f', S.r_overall_mean, S.r_overall_se)
  sprintf('R^2: %.3f ± %.3f', S.R2_mean_mean,  S.R2_mean_se)
  sprintf('EV:  %.3f ± %.3f', S.EV_mean_mean,  S.EV_mean_se)
}, 'FontName','Consolas','VerticalAlignment','top');

disp(S);
