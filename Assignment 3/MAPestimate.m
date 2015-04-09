function [a_MAP, b_MAP, g_obs, u_g, o_g] = MAPestimate(X)

%i.i.d Prior Distribution
u_a = 0;
u_b = 360;
o_a = 10;
o_b = 100;

%Initalization
Parts = X(:,3);
Time = X(:,1) + (X(:,2) - 1)/12;
[m,n] = size(Parts);
%subplot(2,1,1);
%plot(Time,Parts,'b');
%hold on;


%Target function
a_b = sum(Time); 
a_a = sum(Time.^2)+ 1/(o_a^2);
a_const = sum(Time.*Parts); 
b_a = sum(Time);
b_b = m + 1/(o_b^2);
b_const = sum(Parts) + u_b/(o_b^2);

temp = [a_a a_b; b_a b_b]\[a_const;b_const];
a_MAP = temp(1);
b_MAP = temp(2);
%subplot(2,1,1);
%plot(Time,a_MAP*Time+b_MAP,'r')

%Plot error
g_obs = Parts - (a_MAP*Time + b_MAP*ones(m,1));
u_g = sum(g_obs);
o_g = var(g_obs);
%subplot(2,1,2);
%plot(Time,g_obs,'b');








