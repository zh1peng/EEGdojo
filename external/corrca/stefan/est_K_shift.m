function [K, p, W, r, Y, A, gam, r_test] = est_K_shift(data, Nshuffle, gam, data_test)
% estimate number of correlated components through non-parametric statistical test
% based on Nshuffle random circular shifts
%
% Stefan Haufe, 2017

[T D N] = size(data);

if  gam >= 1
  regu = 'tsvd';
else
  regu = 'shrinkage';
end

if isequal(gam, 'auto')
  [~, gam] = shrinkage(reshape(permute(data, [1 3 2]), N*T, [])', []);
end

if T > D
  [W, r, Y, A] = corrca(data, regu, gam);
  if nargin > 3
    [~, r_test] = corrca(data_test, 'fixed', W);
  end
else
  if nargin > 3
    [W, r, Y, ~, r_test] = kcorrca(data, 'linear', 1, data_test, 'mean', min(T, D));
  else
    [W, r, Y] = kcorrca(data, 'linear', 1, [], 'mean', min(T, D));
  end
  A = [];
end

clear r_surro
for ishuffle = 1:Nshuffle
  for in = 1:N
    data_surro(:, :, in) = circshift(data(:, :, in), randi(T));
  end

  if T > D
    [~, r_] = corrca(data_surro, regu, gam);
  else
    [~, r_] = kcorrca(data_surro, 'linear', 1);
  end
  
  r_surro(:, ishuffle) = r_;
end
% toc

p = mean(repmat(r, 1, Nshuffle) < repmat(r_surro(1, :), size(r, 1), 1), 2);

K = sum(p < 0.05);


