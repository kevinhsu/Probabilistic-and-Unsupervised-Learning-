function Kalman_test(X_test, A_TRUE_EM, C_TRUE_EM, Q_TRUE_EM, R_TRUE_EM, A_Ho, C_Ho, Q_Ho, R_Ho, A_Ho_EM, C_Ho_EM, Q_Ho_EM, R_Ho_EM, Acell, Ccell, Qcell, Rcell)

hold on;

X = X_test';

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


%EM at true parameters

[yhat,Vhat,Vjoint,like] = ssm_kalman(X,Y0,Q0,A_TRUE_EM,Q_TRUE_EM,C_TRUE_EM,R_TRUE_EM, 'smooth');
plot(1, sum(like), 'o');

%EM at random parameters
for i = 1:10
    
    A_cell = cell2mat(Acell(i));
    C_cell = cell2mat(Ccell(i));
    Q_cell = cell2mat(Qcell(i));
    R_cell = cell2mat(Rcell(i));
    
    [yhat,Vhat,Vjoint,like] = ssm_kalman(X,Y0,Q0,A_cell,Q_cell,C_cell,R_cell, 'smooth');
    plot(2, sum(like), 'ro');
end

%EM at SSID parameters
[yhat,Vhat,Vjoint,like] = ssm_kalman(X,Y0,Q0,A_Ho_EM,Q_Ho_EM,C_Ho_EM,R_Ho_EM, 'smooth');
plot(3, sum(like), 'g+');

%SSID parameters without EM
[yhat,Vhat,Vjoint,like] = ssm_kalman(X,Y0,Q0,A_Ho,Q_Ho,C_Ho,R_Ho, 'smooth');
plot(4, sum(like), 'y+')