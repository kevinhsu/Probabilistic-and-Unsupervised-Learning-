function [ff] = GPprior()

% Generate some data from a Gaussian Process
X = 1:0.2:20;
X = X';
[M,N] = size(X);
for i = 1:M
    ff(i) = 0 + randn(1,1)*CovKernel(X(i,:),X(i,:));
end

plot(X,ff,'-ob');

%%
function K = CK(N_Time, M_Time)
K = (1+N_Time*M_Time')^3;

function K = CovKernel(N_Time, M_Time)
[N_m,N_n] = size(N_Time);
[M_m,M_n] = size(M_Time);

p1 = 1; %Theta
p2 = 1; %Tao
p3 = 1; %Sigma
p4 = 1; %Pha
p5 = 1; %Nan
p6 = 1; %Cathy
p7 = 1; %Delta

kerf = @(z)((p1^2)*exp(-2*(sin(pi*(z)/p2)).^2/(p3^2)) ...
           +(p1^2)*(p4^2)*exp(-(z).^2/(2*p5^2)) ...
           +(p6^2)*p7);

for i = 1:N_m
    for j = 1:M_m
        s = N_Time(i,:);
        t = M_Time(j,:);       
        K(i,j) = kerf(s-t);
    end
end

K = K + 0.00001*eye(N_m,M_m); %Avoid singular

