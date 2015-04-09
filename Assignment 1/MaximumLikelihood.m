function [p] = MaximumLikelihood(X)
[m, n] = size(X);
p = zeros(n, 1);

for i = 1:n
    sum = 0;
    for j = 1:m
        sum = sum + X(j, i);
    end
    p(i) = sum / m;
end


    