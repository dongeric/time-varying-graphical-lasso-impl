function [best_lambda, best_beta] = tvgl_crossval(train_data, validation_data, rho, norm_type, lambda_set, beta_set)

best_aic = Inf;
for lambda=lambda_set
    for bta = beta_set
        [Thetas, ~] = tvgl_self(train_data, lambda, bta, rho, norm_type);
        
        curr_aic = compute_aic(Thetas, validation_data);
        
        if(curr_aic < best_aic)
            best_lambda = lambda;
            best_beta = bta;
            best_aic = curr_aic;
        end
    end
end
end
                    
                    