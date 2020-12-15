function edge_set = get_random_edges(n)
edge_set = [3; 4];
ind = 1;
for i = 1:n
    for j = 1:n
        if(rand(1) < 0.05 && i ~= j)
           edge_set(:, ind) = [i j];
           ind = ind + 1;
        end
    end
end
end