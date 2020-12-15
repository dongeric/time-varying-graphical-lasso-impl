n = 10;
num_per_epoch = 12;
T = 100;

[cov_series, invcov_series] = generate_cov_matrices(n, T);
train_data = generate_data(cov_series, num_per_epoch);

lambda = 0.2;
beta = 1e3;
rho = 6;
norm_type = 2;

start = tic;
[Thetas, ~] = tvgl_self(train_data, lambda, beta, rho, norm_type);
train_time = toc(start);

[score, precision, recall] = avg_f1(Thetas, invcov_series);