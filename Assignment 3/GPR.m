function [u_Posterior, o_Posterior] = GPR(X)

%Initalization
[a_MAP, b_MAP, g_obs, u_g, o_g] = MAPestimate(X);
Parts = g_obs;
Time = X(:,1) + (X(:,2) - 1)/12;
noise = 1;

N_Time = [];
M_Time = [];

%Generate samples randomly
N_Time = Time;
N_Parts = Parts;
[m,n] = size(N_Time);

M_Time(end+1:end+4,1) = 2007 + ((9:12)' - 1)/12;
for i = 2008:2020
    M_Time(end+1:end+12,1) = i + ((1:12)' - 1)/12;
end
[m2,n2] = size(M_Time);

%Predictions
KMN = CovKernel(M_Time,N_Time);
KNM = CovKernel(N_Time,M_Time);
KNN = CovKernel(N_Time,N_Time);
KMM = CovKernel(M_Time,M_Time);

u_Posterior = KMN * inv(KNN + noise*eye(m)) * (N_Parts);
o_Posterior = KMM - KMN * inv(KNN + noise*eye(m)) * KNM;

for i = 1:m2
    error_SD(i) = sqrt(o_Posterior(i,i) + noise);
end

hold on;
plot(Time,X(:,3),'-b');
plot(M_Time,a_MAP*M_Time+ b_MAP + u_Posterior,'-r');
plot(M_Time,a_MAP*M_Time+ b_MAP + u_Posterior + error_SD', '-g');
plot(M_Time,a_MAP*M_Time+ b_MAP + u_Posterior - error_SD', '-g');
legend('Observed data','Predictive mean','Standard deviation','Standard deviation','SouthEast');

xlabel('Time');
ylabel('Parts per million');

%% sub-functions
function K = CovKernel(N_Time, M_Time)
[N_m,N_n] = size(N_Time);
[M_m,M_n] = size(M_Time);

p1 = 1; %Theta
p2 = 1; %Tao
p3 = 1; %Sigma
p4 = 1; %Pha
p5 = 1; %Nan
p6 = 10; %Cathy
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





