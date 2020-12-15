function aic_val = compute_aic(Thetas, dta)
[n, ~, T] = size(Thetas);

% Compute covariance matrices and number of parameters as done by authors
% of paper, namely looking at which edges change over time
Ss = zeros(n,n,T);
num_params = 0;
for i = 1:T
    if(i == 1)
        num_params = num_params + sum(sum(Thetas(:,:,1) > 0)) ...
            - sum(sum(diag(Thetas(:,:,1)) > 0));
    else
        params_prev = Thetas(:,:,i - 1) > 0;
        params_curr = Thetas(:,:,i) > 0;
        num_params = num_params + sum(sum(bitxor(params_prev, params_curr)));
    end
   
    % Compute covariance matrices
    S = cov(dta(:,:,i));
    Ss(:,:,i) = S;
end


neg_log_likelihood = 0;
for i = 1:T
    neg_log_likelihood = neg_log_likelihood - log_det(Thetas(:,:, i)) ...
        + trace(Ss(:,:,i) * Thetas(:,:,i));
end

aic_val = num_params + neg_log_likelihood / T;
end