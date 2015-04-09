function [p1,p2,p3] = ModelSelection(X)
N = 100;
D = 64;
A = 3;
B = 3;
p1 = 1;
p2 = 1;
p3 = 1;
s = 0;
[m,n] = size(X);

s = ((0.5)^(N)*(((1/beta(A,B))*(0.5)^(A+B-2))))^(1/N); 
p1 = p1 * s; %Model (a)

sum = 0;
for i = 1:m
    for j = 1:n
        if X(i,j) == 1
            sum = sum+1;
        end
    end
end

syms x
s = int((x.^(sum))*((1-x).^(N*D-sum))*((1./beta(A,B))*(x.^(A-1))*((1-x).^(B-1)))^D,0,1);
p2 = double(s^(1/(N*D))); %Model (b)


syms x
for j = 1:n %n=64
    sum = 0;
    for i = 1:m
        if X(i,j) == 1
            sum = sum+1;
        end 
        s = int((x.^(sum))*((1-x).^(N-sum))*((1./beta(A,B))*(x.^(A-1))*((1-x).^(B-1)))^1,0,1);
        p3 = p3 * double(s);
        p3 = p3^(1/(N*D)); %Model (c)
    end
end











