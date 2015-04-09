function [z, u, s, Likelihood] = Binary_EM(K, X, I)

%--------------------------------------------------------------------------
%K: Number of clusters
%X: Observed data
%I: Number of maximum iterations
%z: Responsibility of mixture component
%u: Probability that pixel takes value 1 under mixture component
%s: Mixture proportions
%Likelihood: Array storing log likelihood of P(X|u,s) measured in bits of each step
%--------------------------------------------------------------------------

%Define Variables
[N, D] = size(X);
z = zeros(N,K);
u = zeros(K,D);
s = zeros(K,1);
G = 1/1000000000;
l = 1/G;
l_old = G;

%Random Initialization
for k = 1:K
    for i = 1:D
        u(k,i) = rand(1);
    end
end

s = rand(K,1);
s = s./sum(s);  %Constraint: sum of s should be equal to 1

for iteration = 1:I
     
    %E-step to find responsibility
    for n = 1:N
        for k = 1:K
            mult1 = 1;
            sum2 = 0;
            for i = 1:D
                mult1 = mult1*(u(k,i)^(X(n,i)))*((1 - u(k,i))^(1 - X(n,i)));
            end
            for m = 1:K
                mult21 = 1;
                for i = 1:D
                    mult21 = mult21*(u(m,i)^(X(n,i)))*((1 - u(m,i))^(1 - X(n,i)));
                end
                sum2 = sum2 + s(m)*mult21;
            end
            z(n,k) = (s(k)*mult1)/sum2;        %New z and nomalization
        end
    end
    
    %M-step to find estimated parameters
    for m = 1:K
        sum1 = 0;
        for n = 1:N
            sum1 = sum1 + z(n,m);
        end
        s(m) = sum1./N;          %New s
    end
    
    for m = 1:K
        for j = 1:D
            sum1 = 0;
            sum2 = 0;
            for n =1:N
                sum1 = sum1 + z(n,m)*X(n,j);
            end
            for n = 1:N
                sum2 = sum2 + z(n,m);
            end
            u(m,j) = sum1./sum2;     %New u
        end
    end
    
    %Log Likelihood
    sum2 = 0;
    for n = 1:N
        sum1 = 0;
        for k = 1:K
            mult1 = 1;
            for i = 1:D
                if X(n,i) == 0
                    mult1 = mult1*(1-u(k,i));
                end
                if X(n,i) == 1
                    mult1 = mult1*u(k,i);
                end
            end
            sum1 = sum1 + s(k)*mult1;
        end
        sum2 = sum2 + log(sum1);
    end
    l_old = l;
    l = sum2*log2(exp(1));
    i = iteration;
    Likelihood(i) = l;   %Updated log likelihood measured in bits
    
    %if Convergence, break. Otherwise, iterating until I times
    if abs(l - l_old) < G  
        break;
    end
    
end

plot(Likelihood);











