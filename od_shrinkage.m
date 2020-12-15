function y = od_shrinkage(a, kappa)
[n, ~] = size(a);
a_d = diag(a);
y = max(0, a-kappa) - max(0, -a-kappa);
y(logical(eye(n))) = a_d;
end