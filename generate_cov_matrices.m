function [covariance_series, invcov_series] = generate_cov_matrices(n, T)
% Input: n, the dimension of the data
% Input: T, the number of time steps

% Generate one inverse covariance
inv_cov = zeros(n);
while(det(inv_cov) < 1e-3 || ~all(diag(inv(inv_cov)) > 0) )
    S = sprand(n, n, 0.2);

    inv_cov = full(S * S.');
    
end

[U, S, V] = svd(inv_cov);

for i = 1:3
    S(i,i) = S(i,i) * 10;
end

inv_cov = U * S * V';

invcov_series = zeros(n, n, T);
covariance_series = zeros(n, n, T);

cov = inv(inv_cov);
for i = 1:T/2
    invcov_series(:, :, i) = inv_cov;
    covariance_series(:, :, i) = cov;
end
    
% Generate second inverse covariance
inv_cov = zeros(n);
while(det(inv_cov) < 1e-3 || ~all(diag(inv(inv_cov)) > 0) )
    S = sprand(n, n, 0.2);

    inv_cov = full(S * S.');
end

[U, S, V] = svd(inv_cov);

for i = 1:3
    S(i,i) = S(i,i) * 10;
end

inv_cov = U * S * V';

cov = inv(inv_cov);
for i = T/2+1:T
    invcov_series(:, :, i) = inv_cov;
    covariance_series(:, :, i) = cov;
end

invcov_series = invcov_series .* (invcov_series > 1e-3);

end

