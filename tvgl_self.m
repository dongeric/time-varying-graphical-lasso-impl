function [Thetas, history] = tvgl_self(D, lambda, bta, rho, norm_type)
% code skeleton taken from https://web.stanford.edu/~boyd/papers/admm/
QUIET    = 0;
MAX_ITER = 5e2;
ABSTOL   = 1e-3;
RELTOL   = 1e-4;

[num_per_epoch, n, T] = size(D);
Thetas = zeros(n, n, T);
best_Thetas = zeros(n, n, T);
curr_best_val = Inf;

Z0s = zeros(n, n, T);
Z1s = zeros(n, n, T - 1);
Z2s = zeros(n, n, T - 1);

U0s = zeros(n, n, T);
U1s = zeros(n, n, T - 1);
U2s = zeros(n, n, T - 1);


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


if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end


for k = 1:MAX_ITER
    % theta-update
    for i = 1:T
        Si = Ss(:,:,i);
        if(i > 1 && i < T)
            A = (Z0s(:,:,i) - U0s(:,:,i) + Z1s(:,:,i) - U1s(:,:,i) + ...
                 Z2s(:,:,i - 1) - U2s(:,:,i- 1)) / 3;
            A = (A + A.') / 2;
            eta = num_per_epoch / (3 * rho); 
        elseif (i == 1)
            A = (Z0s(:,:,i) - U0s(:,:,i) + Z1s(:,:,i) - U1s(:,:,i)) / 2;
            A = (A + A.') / 2;
            eta = num_per_epoch / (2 * rho);
        elseif (i == T)
            A = (Z0s(:,:,i) - U0s(:,:,i) + Z2s(:,:,i - 1) - U2s(:,:,i - 1)) / 2;
            A = (A + A.') / 2;
            eta = num_per_epoch / (2 * rho);
        end
        [Q,L] = eig(A / eta - Si);
        es = diag(L);
        xi = (es + sqrt(es.^2 + 4 / eta));
        Thetas(:,:,i) = Q*diag(xi)*Q' * eta / 2;
    end

    % z0-update
    Z0olds = Z0s;
    for i = 1:T
        Z0s(:,:,i) = od_shrinkage(Thetas(:,:,i) + U0s(:,:,i), lambda/rho);
    end
    
    %z1 and z2 - update
    Z1olds = Z1s;
    Z2olds = Z2s;
    for i = 2:T
       A = (Thetas(:,:,i - 1) - Thetas(:,:,i) + U1s(:,:, i - 1) - U2s(:,:,i - 1));
       eta = 2 * bta / rho;
       
       if(norm_type == 1) % L1 norm
           E = shrinkage(A, eta);
       elseif(norm_type == 2) % L2 norm
           E = block_shrinkage(A, eta);
       elseif(norm_type == 3) % Laplacian Norm
           E = A ./ (1 + 2 * eta);
       elseif(norm_type == 6) % New: Nuclear Norm
           [U, S, V] = svd(A);
           E = U * shrinkage(S, eta) * V';
       end
       
       core = Thetas(:,:,i) + Thetas(:,:,i - 1) + U1s(:,:, i - 1) + U2s(:,:,i - 1);
       Z1s(:,:,i - 1) = (core - E) / 2;
       Z2s(:,:,i - 1) = (core + E) / 2;
    end
    

    % Dual updates
    U0s = U0s + Thetas - Z0s;
    U1s = U1s + Thetas(:,:,1:T-1) - Z1s;
    U2s = U2s + Thetas(:,:,2:T) - Z2s;
    
    
    % diagnostics, reporting, termination checks
    
    % primal norm
    r = sum(sum(sum((Thetas - Z0s).^2))) + ...
        sum(sum(sum((Thetas(:,:,1:T-1) - Z1s).^2))) + ...
        sum(sum(sum((Thetas(:,:,2:T) - Z2s).^2)));
    
    % dual norm
    s = sum(sum(sum((Z0olds - Z0s).^2))) + ...
        sum(sum(sum((Z1olds - Z1s).^2))) + ...
        sum(sum(sum((Z2olds - Z2s).^2)));
    
    

    history.objval(k)  = tv_objective(Ss, Thetas, lambda, bta, norm_type);
    
    if(history.objval(k) < curr_best_val)
        curr_best_val = history.objval(k);
        best_Thetas = Thetas;
    end

    history.r_norm(k)  = sqrt(r);
    history.s_norm(k)  = sqrt(s);

    history.eps_pri(k) = sqrt(n*n*T)*ABSTOL + ...
        RELTOL*max(sqrt(sum(sum(sum(Thetas.^2)))), sqrt(sum(sum(sum(Z0s.^2)))));
    
    history.eps_dual(k)= sqrt(n*n*T)*ABSTOL + ...
        RELTOL*sqrt(rho * (sum(sum(sum(U0s.^2)))));


    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end
end

Thetas = best_Thetas .* (best_Thetas > 1e-3);
end