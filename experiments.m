% Get Synthetic Data (Train, Validation, Test)

n = 10;
num_per_epoch = 12;
T = 100;

[cov_series, invcov_series] = generate_cov_matrices(n, T);

train_data = generate_data(cov_series, num_per_epoch);

validation_data = generate_data(cov_series, num_per_epoch);

test_data = generate_data(cov_series, num_per_epoch);

% Cross Validate to get parameters for each norm

lambda_set = [0.04, 0.16, 0.64];

beta_set = [1e3, 1e4, 1e6];

rho_1 = 1;
% l1 norm parameters
[lambda_1, beta_1] = tvgl_crossval(train_data, validation_data, ...
                        rho_1, 1, lambda_set, beta_set)
                    

rho_2 = 1;
% l2 norm parameters
[lambda_2, beta_2] = tvgl_crossval(train_data, validation_data, ...
                        rho_2, 2, lambda_set, beta_set)
                    

rho_3 = 1;
% laplacian norm parameters
[lambda_3, beta_3] = tvgl_crossval(train_data, validation_data, ...
                        rho_3, 3, lambda_set, beta_set)

rho_6 = 1;
% nuclear norm parameters
[lambda_6, beta_6] = tvgl_crossval(train_data, validation_data, ...
                        rho_6, 6, lambda_set, beta_set)

% Run tests for each norm

start1 = tic;
[Thetas_1, ~] = tvgl_self(test_data, lambda_1, beta_1, rho_1, 1);
l1_train_time = toc(start1);

[score_1, precision_1, recall_1] = avg_f1(Thetas_1, invcov_series);

start2 = tic;
[Thetas_2, ~] = tvgl_self(test_data, lambda_2, beta_2, rho_2, 2);
l2_train_time = toc(start2);

[score_2, precision_2, recall_2] = avg_f1(Thetas_2, invcov_series);

start3 = tic;
[Thetas_3, ~] = tvgl_self(test_data, lambda_3, beta_3, rho_3, 3);
l3_train_time = toc(start3);

[score_3, precision_3, recall_3] = avg_f1(Thetas_3, invcov_series);

start6 = tic;
[Thetas_6, ~] = tvgl_self(test_data, lambda_6, beta_6, rho_6, 6);
l6_train_time = toc(start6);

[score_6, precision_6, recall_6] = avg_f1(Thetas_6, invcov_series);


