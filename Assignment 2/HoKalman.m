function [A,Q,C,R,Likelihood] = HoKalman(X, L)

%L is maximum lag

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

[D,T] = size(X); %5*1000
[K,K] = size(A); %4*1000

maxlag = L;

% Approxiamtion approach 1
% M = cell(2*maxlag - 1,1);
% for tt = 1:2*maxlag - 1
%     sum1 = 0;
%     for t = maxlag + 1:T - 2*maxlag + 1
%         sum1 = sum1 + X(:,tt+t)*X(:,t)';
%     end
%     M(tt) = {(1/T)*sum1};
% end
% H = cell2mat(M(hankel([1:maxlag], [maxlag:2*maxlag-1])));

%Approxiamtion approach 2
H = zeros(D*maxlag, D*maxlag);
M = zeros(D*maxlag, D*maxlag);
for t = maxlag+1:T-maxlag+1
    X_positive = [];
    X_negative = [];
    for m = 0:maxlag-1
        X_positive = [X_positive,X(:,t+m)'];
        X_negative = [X_negative,X(:,t-m-1)'];
    end
    M = M + X_positive'*X_negative;
end
H = (1/(T-2*maxlag + 1))*M;

[Xi,SV,Ups] = svds(H,K);
Xi = Xi*sqrt(SV);
Ups = sqrt(SV)*Ups';

Ahat = Xi(1:end-D,:)\Xi(D+1:end,:);
Chat = Xi(1:D,:);
PIhat = Ahat\(Ups(:,1:D)/Chat');
Qhat = PIhat - Ahat*PIhat*Ahat';

sum1 = 0;
for t = 1:T
    sum1 = sum1 + X(:,t)*X(:,t)';
end
Rhat = (1/T)*sum1 - Chat*PIhat*Chat';

A = Ahat;
Q = Qhat;
C = Chat;
R = Rhat;

[yhat,Vhat,Vjoint,like] = ssm_kalman(X,Y0,Q0,A,Q,C,R, 'smooth');

Likelihood = sum(like);






