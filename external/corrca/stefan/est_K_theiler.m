function K = est_K_theiler(data, Nshuffle, gam)
% estimate number of correlated components through non-parametric statistical test
% based on Nshuffle random phase scramblings 
% (Theiler et al. 1992, Physica D; Pritchard and Theiler 1994, PRL)
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
  [~, r] = corrca(data, regu, gam);
else
  [~, r] = kcorrca(data, 'linear', 1);
end


data_gauss = sort(repmat(std(data), T, 1, 1).*randn(size(data))+repmat(mean(data), T, 1, 1));
[so1 ing1] = sort(data);
[so2 ing2] = sort(ing1);
for in = 1:N
  for id = 1:D
    data_gauss(:, id, in) = data_gauss(ing2(:, id), id, in);
  end
end
data_gauss_f = fft(data_gauss);
[phase, amp] = cart2pol(real(data_gauss_f), imag(data_gauss_f));


clear r_surro
for ishuffle = 1:Nshuffle
%   ishuffle
  
  phase_ = repmat(2*pi*rand(T/2-1, 1, N)-pi, 1, D, 1);
  phase(2:T/2, :, :) = phase_;
  phase((T/2+2):end, :, :) = -flipud(phase_);
  [rea, ima] = pol2cart(phase, amp);
  data_gauss_f_surro = rea + sqrt(-1)*ima;
  data_gauss_surro = real(ifft(data_gauss_f_surro));

  [so3 ing3] = sort(data_gauss_surro);
  [so4 ing4] = sort(ing3);
  
  data_surro = [];
  for in = 1:N
    for id = 1:D
      data_surro(:, id, in) = so1(ing4(:, id, in), id, in);
    end
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


