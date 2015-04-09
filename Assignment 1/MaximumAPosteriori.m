function [p] = MaximumAPosteriori(X)

Alpha = 3;
Beta = 3;

[m, n] = size(X);

for i = 1:n
    sum = 0;
    for j = 1:m
        sum = sum + X(j, i);
    end
    p(i) = (sum + (Alpha - 1)) / (m + Alpha + Beta - 2);
end

        