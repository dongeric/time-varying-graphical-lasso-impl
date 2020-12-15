function tv_data = generate_data(cov_series, num_per_epoch)

[n, ~, T] = size(cov_series);

tv_data = zeros(num_per_epoch, n, T);
mu = zeros(n, 1);

for i = 1:T
    cov = cov_series(:,:,i);
    
    x = mvnrnd(mu, cov, num_per_epoch);
    
    tv_data(:, :, i) = x;
end

end