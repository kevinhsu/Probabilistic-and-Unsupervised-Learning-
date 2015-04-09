function [A,Q,C,R,Likelihood] = HoKalman_EM(X)

%Initialization
X = X'; %5*1000
Y0 = zeros(4,1);
Q0 = eye(4,4);

A = [0.9627   -0.0424    0.0396   -0.0279;
    0.0495    0.9818   -0.0503   -0.0246;
   -0.0404    0.0103    1.0009   -0.0093;
    0.0244    0.0397    0.0057    0.9993];

Q = [  0.0258   -0.0087    0.0060    0.0003;
   -0.0096    0.0374    0.0158   -0.0056;
    0.0072    0.0158    0.0197   -0.0020;
    0.0029   -0.0055   -0.0032    0.0187];
C = [ -1.0919   -1.0073   -0.0885   -0.3143;
   -0.8138    0.6903    0.4026   -0.3957;
   -1.0250    0.1831   -0.4575    0.0527;
   -1.1346   -0.0557    0.5601    0.4311;
   -0.9265   -0.1341    0.2390   -0.4169];
R = [
    0.9553    0.0051    0.0202    0.0431    0.0469;
    0.0029    1.0197    0.0189    0.0049   -0.0318;
    0.0106   -0.0282    1.0889    0.0001    0.0703;
    0.0723   -0.0116    0.0008    1.0348    0.0292;
    0.0462   -0.0249    0.1037    0.0238    0.9945];

[S,T] = size(X);
M = 50;
Likelihood = zeros(M,0);
[yhat,Vhat,Vjoint,like] = ssm_kalman(X,Y0,Q0,A,Q,C,R, 'smooth');

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
    Likelihood(i,1) = sum(like);
       
end

C = C_new;
A = A_new;
R = R_new;
Q = Q_new;

plot(Likelihood);