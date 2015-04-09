function [A, Q, C, R, Likelihood] = Kalman_EM(X)

%Initialization
X = X'; %5*1000
Y0 = zeros(4,1);
Q0 = eye(4,4);
A = 0.99*[cos(2*pi/180) -sin(2*pi/180) 0 0;
    sin(2*pi/180) cos(2*pi/180) 0 0;
    0 0 cos(2*pi/180) -sin(2*pi/90);
    0 0 sin(2*pi/90) cos(2*pi/90)];
Q = eye(4,4) - A*A';
C = [1 0 1 0;
    0 1 0 1;
    1 0 0 1;
    0 0 1 1;
    0.5 0.5 0.5 0.5];
R = eye(5,5);

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

plot(Likelihood);












