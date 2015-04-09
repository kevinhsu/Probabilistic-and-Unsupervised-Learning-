function [prior,prediction] = Gaussian_Process( kernel,xi2,X,Y,X_pre )
% Gaussian_Process: providing prior,posterior and prediction with given
% training data and kernel according to gaussian process
% input:
% kernel:     a string of kernel name
% xi2:        observation noise
% X:          training input
% Y:          observed value
% X_pre:      prediction input
% output:
% prior:      prior of random sample
% prediction: prediction mean(1st row), standard deviation(2nd row)

% computing prediction
if nargin>4
    K=feval(kernel,X,X,xi2);
    K=(K+K')/2;
    prior=[];
    K_pre_X=feval(kernel,X_pre,X);
    ave=K_pre_X/K*Y';
    std_err=sqrt(diag(feval(kernel,X_pre,X_pre,xi2)-K_pre_X/K*K_pre_X'));
    prediction=[ave';std_err'];
% computing prior
else
    [~,n]=size(X);
    prediction=[];
    K=feval(kernel,X,X);
    K=(K+K'/2);
    L=chol(K+0.0001*eye(n));                    % avoid numerical problem
    prior=randn(3,n)*L;                         % generate 3 different sample prior functions
end   
end
