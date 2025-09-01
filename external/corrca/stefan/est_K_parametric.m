function K = est_K_parametric(data, Nsplits, gam)
% estimate number of correlated components through parametric statistical test. 
% This is done on test data for Nsplits splits of the data into training and 
% test sets. On the test sets, an F distribution is assumed for the ISC
% (see corrca.c). The number of significant components is determined using
% Bonferroni correction.
%
% Stefan Haufe, 2017

[T D N] = size(data);

if  gam >= 1
  regu = 'tsvd';
else
  regu = 'shrinkage';
end

clear K
for isplit = 1:Nsplits
  
  pe = randperm(T);
  
  data_train = data(pe(1:ceil(T/2)), :, :);
  data_test = data(pe((ceil(T/2)+1):end), :, :);
  
  if isequal(gam, 'auto')
    [~, gam] = shrinkage(reshape(permute(data_train, [1 3 2]), N*ceil(T/2), [])', []);
  end
  
  if T > D
    W_train = corrca(data_train, regu, gam);
    [~, ~, ~, ~, p] = corrca(data_test, 'fixed', W_train);
  else
    [W,ISC,Y,Ytest,r] = kcorrca(data_train, 'linear', 1, data_test, 'mean', min(floor(T/2), D));
  end
 
  K(isplit) = sum(p < 0.05/D);

end

K = median(K);


