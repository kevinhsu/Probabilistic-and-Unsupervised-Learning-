function [Acell,Ccell,Qcell,Rcell,Likelihoodcell] = Kalman_EM_10(X)

X = X';
Y0 = zeros(4,1);
Q0 = eye(4,4);

[D, T] = size(X);

K = D - 1;

Acell = cell(10,1);
Ccell = cell(10,1);
Qcell = cell(10,1);
Rcell = cell(10,1);
Likelihoodcell = cell(10,1);

for trial = 1:10
    
    S = diag(rand(K,1));
    U = orth(rand(K,K));
    A = U' * S * U;
    
    S = diag(rand(K,1));
    U = orth(rand(K,K));
    Q = U' * S * U;
    Q = eye(K,K) - Q*Q';   

    S = diag(rand(D,1));
    U = orth(rand(D,D));
    R = U' * S * U;

    C = rand(D,K);
    for i = 1:D;
        C(i,:) = C(i,:)/sum(C(i,:));
    end
    
    [S,T] = size(X);
    M = 50;
    Likelihood = zeros(M+1,0);
    [yhat,Vhat,Vjoint,like] = ssm_kalman(X,Y0,Q0,A,Q,C,R, 'smooth');
    Likelihood(1,1) = sum(like);
    
    for iteration = 1:M
        
        %C_new
        sum1 = 0;
        sum2 = 0;
        for t = 1:T
            sum1 = sum1 + X(:,t)*yhat(:,t)';
            sum2 = sum2 + yhat(:,t)*yhat(:,t)' + cell2mat(Vhat(t));
        end
        C_new = sum1*(sum2)^(-1);
        %A_new
        sum1 = 0;
        sum2 = 0;
        for t = 2:T
            sum1 = sum1 + yhat(:,t)*yhat(:,t-1)' + cell2mat(Vjoint(t-1));
            sum2 = sum2 + yhat(:,t-1)*yhat(:,t-1)' + cell2mat(Vhat(t-1));
        end
        A_new = sum1*(sum2)^(-1);
        %R_new
        sum1 = 0;
        sum2 = 0;
        for t = 1:T
            sum1 = sum1 + X(:,t)*X(:,t)';
            sum2 = sum2 + X(:,t)*yhat(:,t)'*C_new';
        end
        R_new = (1/T)*(sum1 - sum2);
        %Q_new
        sum1 = 0;
        sum2 = 0;
        for t = 2:T
            sum1 = sum1 + yhat(:,t)*yhat(:,t)' + cell2mat(Vhat(t));
            sum2 = sum2 + (yhat(:,t)*yhat(:,t-1)' + cell2mat(Vjoint(t-1)))*A_new';
        end
        Q_new = (1/(T-1))*(sum1 - sum2);
        
        [yhat, Vhat, Vjoint, like] = ssm_kalman(X,Y0,Q0,A_new,Q_new,C_new,R_new, 'smooth');
        i = iteration;
        Likelihood(i+1,1) = sum(like);
        
    end
    
    C = C_new;
    A = A_new;
    R = R_new;
    Q = Q_new;
    
    i = trial;
    
    Ccell(i) = {C};
    Acell(i) = {A};
    Rcell(i) = {R};
    Qcell(i) = {Q};
    Likelihoodcell(i) = {Likelihood}; 
       
end

for trial = 1:10
   i = trial;
   subplot(5,2,i);
   plot(cell2mat(Likelihoodcell(i))); 
end


