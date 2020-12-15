function y = block_shrinkage(a, kappa)
    thresh = max(vecnorm(a)-kappa, 0);
    y = thresh.*(1-kappa./vecnorm(a)).*a;
end