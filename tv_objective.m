function obj = tv_objective(Ss, Thetas, lambda, bta, norm_type)
[n, ~, T] = size(Thetas);


Z0s = Thetas;
Z1s = Thetas(:,:,1:T-1);
Z2s = Thetas(:,:,2:T);
obj = 0;
for i = 1:T
    obj = obj - log_det(Thetas(:,:, i)) + trace(Ss(:,:,i) * Thetas(:,:,i)) ...
        + lambda * (sum(sum(abs(Z0s(:,:,i)))) - trace(abs(Z0s(:,:,i))));
end

if(norm_type == 1) % l1 norm
    for j = 2:T
        obj = obj + bta * sum(sum(sum(abs(Z1s(:,:,j-1) - Z2s(:,:,j-1)))));
    end 
elseif(norm_type == 2) % l2 norm
    for j = 2:T
        obj = obj + bta * sum(sum(vecnorm(Z1s(:,:,j-1) - Z2s(:,:,j-1))));
    end
elseif(norm_type == 3) % Laplacian
    for j = 2:T
        obj = obj + bta * sum(sum(sum((Z1s(:,:,j-1) - Z2s(:,:,j-1)).^2)));
    end 
elseif(norm_type == 6) % New: nuclear norm
    for j = 2:T
        obj = obj + bta * norm_nuc(Z1s(:,:,j-1) - Z2s(:,:,j-1));
    end 
end

end