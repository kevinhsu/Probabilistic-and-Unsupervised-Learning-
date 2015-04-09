function parameters(X)

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

logdet = @(A)(2*sum(log(diag(chol(A)))));
[Y,V,~,L] = ssm_kalman(X',Y0,Q0,A,Q,C,R, 'filt');
subplot(2,2,1);
plot(Y');
subplot(2,2,2);
plot(cellfun(logdet,V));

[Y,V,Vj,L] = ssm_kalman(X',Y0,Q0,A,Q,C,R, 'smooth');
subplot(2,2,3);
plot(Y');
subplot(2,2,4);
plot(cellfun(logdet,V));