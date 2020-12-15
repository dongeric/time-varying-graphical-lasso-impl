function matrices = tvgl_cvx(D, lambda, bta, norm_type)
tv_data = D;
[num_per_epoch, n, T] = size(tv_data)

% Compute covariance matrices
Ss = zeros(n,n,T);
for i = 1:T 
    if(num_per_epoch == 1)
        x = D(:,:,i);
        S = x.' * x;
    else
        S = cov(D(:,:,i));
    end
    Ss(:,:,i) = S;
end

cvx_begin quiet sdp
    obj = 0;
    variable thetas(n, n, T)
    for i = 1:T
        obj = obj - log_det(thetas(:,:, i)) + trace(Ss(:,:,i) * thetas(:,:, i)) ...
            + lambda * (sum(sum(abs(thetas(:,:,i)))));
    end
    if(norm_type == 1) % l1 norm
        for j = 2:T
            obj = obj + bta * sum(sum(sum(abs(thetas(:,:,j) - thetas(:,:,j-1)))));
        end 
    elseif(norm_type == 2) % l2 norm
        for j = 2:T
            obj = obj + bta * sum(sum(norm(thetas(:,:,j) - thetas(:,:,j-1), 2)));
        end
    elseif(norm_type == 3) % Laplacian
        for j = 2:T
            obj = obj + bta * sum(sum(sum((thetas(:,:,j) - thetas(:,:,j-1)).^2)));
        end 
    elseif(norm_type == 6) % New: nuclear norm
        for j = 2:T
            obj = obj + bta * norm_nuc((thetas(:,:,j) - thetas(:,:,j-1)));
        end 
    end
    minimize obj
    subject to
    for i = 1:T
        thetas(:,:,i) >= 0;
    end

cvx_end

matrices = thetas .* (thetas > 1e-7);
end