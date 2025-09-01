function run_sim8
% Gaussian white signal and noise time series,
% spatial correlation structure of signal/noise either the same or 
% different across subjects
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
T = 200;

% K: number of correlated components
K = 10;

% SNR
db_snr = 0;
snr = 10.^(db_snr/20)./(1+10.^(db_snr/20)); %=0.5
% snr = 0.1;
% dbs = 20*log10(snrs./(1-snrs));

% regularization parameter
gam = 0;

distrs = {'same_signal_same_noise', 'same_signal_different_noise', ...
  'different_signal_same_noise', 'different_signal_different_noise'};

clear K_shift K_phase K_param K_spectral Y_acc ISC_av ISC_av_test ISC_test Y_subspace Y_subspace_regr
for idistr = 1:length(distrs)

  for irep = 1:Nrep
    
    switch idistr
      case 1
        [data, Y_true, A_true] = gen_data(2*T, D, N, ones(1, K), snr, 0, 1, 1, 0, 'Gauss');
      case 2
        [data, Y_true, A_true] = gen_data(2*T, D, N, ones(1, K), snr, 0, 1, 0, 0, 'Gauss');
      case 3
        [data, Y_true, A_true] = gen_data(2*T, D, N, ones(1, K), snr, 0, 0, 1, 0, 'Gauss');
      case 4
        [data, Y_true, A_true] = gen_data(2*T, D, N, ones(1, K), snr, 0, 0, 0, 0, 'Gauss');
    end
    
    data_test = data(T+1:end, :, :);
    data = data(1:T, :, :);
    Y_true = Y_true(1:T, :, :);
    
    tic
    [K_shift(idistr, irep), p, W, ISC(:, idistr, irep), Y, A, gam, ISC_test(:, idistr, irep)] = ...
      est_K_shift(data, Nshuffle, gam, data_test);
    t1 = toc;

    tic
    K_phase(idistr, irep) = est_K_theiler(data, Nshuffle, gam);
    t2 = toc;
    
    tic
    K_param(idistr, irep) = est_K_parametric(data, Nsplits, gam);
    t3 = toc;
    
    K_spectral(idistr, irep) = sum(sign(ISC(:, idistr, irep)));

%     Y_est = reshape(permute(Y(:, 1:K, :), [1 3 2]), T*N, K);
%     Y_true = reshape(permute(Y_true, [1 3 2]), T*N, K);
    Y_est = mean(Y(:, 1:K, :), 3);
    Y_true = mean(Y_true, 3);
    Y_rot = Y_est*(Y_est\Y_true);
    Y_acc(idistr, irep) = mean(max(abs(corr(Y_true, Y_est))'));
    Y_subspace(idistr, irep) = 1-subspace(Y_est, Y_true)/(pi/2);
    Y_subspace_regr(idistr, irep) = mean(diag(corr(Y_rot, Y_true)));
    
    ISC_av(idistr, irep) = mean(ISC(1:K, idistr, irep));
    ISC_av_test(idistr, irep) = mean(ISC_test(1:K, idistr, irep)); 

    [idistr irep t1 t2 t3]
  end
end

save('mat/sim8a_results', 'K_shift', 'K_phase', 'K_param', 'K_spectral', 'Y_acc', ...
  'Y_subspace', 'Y_subspace_regr', 'ISC_av', 'ISC_av_test', ...
  'ISC', 'ISC_test', 'Nrep', 'Nshuffle', 'Nsplits', 'N', 'D', 'T', 'K', 'db_snr', 'snr', 'distrs', 'gam')
  
