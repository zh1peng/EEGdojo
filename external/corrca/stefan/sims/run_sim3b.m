function run_sim3b
% colored (pink) Gaussian signal and noise
% vary T
%
% Stefan Haufe, 2017

% Nrep: number of repetitions to obtain confidence interval
Nrep = 100;

% Nshuffle: number of phase shuffles/circular shifts for non-parametric
% significance test
Nshuffle = 1000;

% number of train/test splits for parametric significance test
Nsplits = 100;

% N: number of subjects/viewings
N = 5;

% D: number of data channels
D = 30;

% T: number of time samples
Ts = [100 300 1000 3000 10000];

% K: number of correlated components
K = 10;

% SNR
db_snr = 0;
snr = 10.^(db_snr/20)./(1+10.^(db_snr/20)); %=0.5
% snr = 0.1;
% dbs = 20*log10(snrs./(1-snrs));

% regularization parameter
gam = 0;


clear K_shift K_phase K_param K_spectral A_acc Y_acc ISC_av ISC_av_test ISC_test A_subspace A_subspace_regr Y_subspace Y_subspace_regr
for it = 1:length(Ts)
  T = Ts(it);
  
  for irep = 1:Nrep
    
    [data, Y_true, A_true] = gen_data(2*T, D, N, ones(1, K), snr, 0, 1, 1, 1, 'Gauss');

    data_test = data(T+1:end, :, :);
    data = data(1:T, :, :);
    Y_true = Y_true(1:T, :, :);
    
    tic
    [K_shift(it, irep), p, W, ISC(:, it, irep), Y, A, gam, ISC_test(:, it, irep)] = ...
      est_K_shift(data, Nshuffle, gam, data_test);
    t1 = toc;
    
    tic
    K_phase(it, irep) = est_K_theiler(data, Nshuffle, gam);
    t2 = toc;
    
    tic
    K_param(it, irep) = est_K_parametric(data, Nsplits, gam);
    t3 = toc;
    
    K_spectral(it, irep) = sum(sign(ISC(:, it, irep)));

    A_acc(it, irep) = mean(max(abs(corr(A(:, 1:K), A_true'))));
    A_subspace(it, irep) = 1-subspace(A(:, 1:K), A_true')/(pi/2);
    A_rot = A(:, 1:K)*(A(:, 1:K)\A_true');
    A_subspace_regr(it, irep) = mean(diag(corr(A_rot, A_true')));
    
%     Y_est = reshape(permute(Y(:, 1:K, :), [1 3 2]), T*N, K);
%     Y_true = reshape(permute(Y_true, [1 3 2]), T*N, K);
    Y_est = mean(Y(:, 1:K, :), 3);
    Y_true = mean(Y_true, 3);
    Y_rot = Y_est*(Y_est\Y_true);
    Y_acc(it, irep) = mean(max(abs(corr(Y_true, Y_est))'));
    Y_subspace(it, irep) = 1-subspace(Y_est, Y_true)/(pi/2);
    Y_subspace_regr(it, irep) = mean(diag(corr(Y_rot, Y_true)));
  
    ISC_av(it, irep) = mean(ISC(1:K, it, irep));
    ISC_av_test(it, irep) = mean(ISC_test(1:K, it, irep));
    
    [it irep t1 t2 t3]
  end
end

save('mat/sim3b_results', 'K_shift', 'K_phase', 'K_param', 'K_spectral', 'A_acc', 'A_subspace', 'A_subspace_regr', ...
  'Y_subspace', 'Y_subspace_regr', 'Y_acc', 'ISC_av', 'ISC_av_test', ...
  'ISC', 'ISC_test', 'Nrep', 'Nshuffle', 'Nsplits', 'N', 'D', 'Ts', 'K', 'db_snr', 'snr', 'gam')
  
