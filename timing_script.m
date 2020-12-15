lambda = 0.4;
beta = 1e3;
rho = 6;

n_set = [6 8 10];
T_set = [20 50 100];
norm_types = [1 2 3 6];

cvx_results = zeros(length(n_set), length(norm_types));
tvgl_results = zeros(length(n_set), length(norm_types));

for i = 1:length(n_set)
    n = n_set(i)
    T = T_set(i)
    num_per_epoch = 12;

    [cov_series, invcov_series] = generate_cov_matrices(n, T);
    train_data = generate_data(cov_series, num_per_epoch);
    
    for j = 1:length(norm_types)
         norm_type = norm_types(j)
        start = tic;
        [Thetas, ~] = tvgl_self(train_data, lambda, beta, rho, norm_type);
        train_time = toc(start)
        
        tvgl_results(i, j) = train_time;
        
        start = tic;
        tvgl_cvx(train_data, lambda, beta, norm_type);
        train_time = toc(start)

        cvx_results(i, j) = train_time;
    end
end

%% Last part only for tvgl

tvgl_high_results = zeros(1, length(norm_types));

num_per_epoch = 12;
n = 15;
T = 200;

[cov_series, invcov_series] = generate_cov_matrices(n, T);
train_data = generate_data(cov_series, num_per_epoch);

for j = 1:length(norm_types)
        norm_type = norm_types(j)
        start = tic;
        [Thetas, ~] = tvgl_self(train_data, lambda, beta, rho, norm_type);
        train_time = toc(start)
        
        tvgl_high_results(1, j) = train_time;
end
