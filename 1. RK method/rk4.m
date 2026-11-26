% Use the RK (direct) method for the system to approximate the solution of
% the following second order differential equation and compare the result
% to the actual solution.
%
% t^2*y'' - 2*t*y' + 2*y = t^3*ln(t),   1 <= t <= 2,
% y(1) = 1,  y'(1) = 0,  with h = 0.2
%
% Actual solution:
% y(t) = (7/4)*t + (1/2)*t^3*ln(t) - (3/4)*t^3

function rk4
clc
clf
clear all
a=1;
b=2;
h=0.2;
n=(b-a)/h;
t=linspace(a,b,n+1);
alpha1=1;
alpha2=0;
u(1)=alpha1;
u1(1)=alpha2;
for i=1:n
    k1=h*f(t(i),u(i),u1(i));
    k2=h*f(t(i)+h/2,u(i)+u1(i)*h/2+h*k1/8,u1(i)+k1/2);
    k3=h*f(t(i)+h/2,u(i)+u1(i)*h/2+h*k1/8,u1(i)+k2/2);
    k4=h*f(t(i)+h,u(i)+u1(i)*h+h*k3/2,u1(i)+k3);
    u(i+1)=u(i)+h*(u1(i)+(k1+k2+k3)/6);
    u1(i+1)=u1(i)+(k1+2*k2+2*k3+k4)/6;
end
plot(t,u,'-b');
hold on
eu=7/4*t+1/2*t.^3.*log(t)-3/4*t.^3;
stem(t,eu,'*k');
hold on
legend('approximate solution','actual solution');
fprintf(' t approximate exact error \n');
for i=1:n+1
    fprintf(' %11.5f %11.5f %11.5f %11.5f \n',t(i),u(i),eu(i),abs(eu(i)-u(i)));
end
function z=f(t,u,u1)
z=2/t*u1-2*u/t^2+t*log(t);