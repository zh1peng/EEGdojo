%% no shift
x = randn(3000, 1, 1);
[b a] = butter(5, 0.3);
x = filter(b, a, x);

z = repmat(x, 1, 10, 20);
z = z ./ norm(z(:));

noise = randn(size(z));
noise = noise ./ norm(noise(:));

meas = z + noise;

[W,ISC,A,Y] = corrca(meas);
ISC


%% repeats shifted
shft = randi(21, 1, 1, 20)-10;
for jj = 1:20
  z(:, :, jj) = repmat(circshift(x, shft(1, 1, jj), 1), 1, 10, 1);
end

meas = z + noise;

[W,ISC,A,Y] = corrca(meas);
ISC

%% repeats and channels shifted
shft = randi(21, 1, 10, 20)-10;
for ii = 1:10
  for jj = 1:20
    z(:, ii, jj) = circshift(x, shft(1, ii, jj), 1);
  end
end

meas = z + noise;

[W,ISC,A,Y] = corrca(meas);
ISC

%% no shift, but half of the repeats has a sign flip
z = repmat(x, 1, 10, 20);
z(:, :, 11:20) = -z(:, :, 11:20);
z = z ./ norm(z(:));

noise = randn(size(z));
noise = noise ./ norm(noise(:));

meas = z + noise;

[W,ISC,A,Y] = corrca(meas);
ISC