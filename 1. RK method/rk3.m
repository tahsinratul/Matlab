% Use the RK method for the system to approximate the solution of the
% following second order differential equation and compare the result
% to the actual solution.
%
% t^2*y'' - 2*t*y' + 2*y = t^3*ln(t),   1 <= t <= 2,   h = 0.1,
% y(1) = 1,  y'(1) = 0
%
% Actual solution:
% y(t) = (7/4)*t + (1/2)*t^3*ln(t) - (3/4)*t^3

function rk3
clc
clf
clear all
a=1;
b=2;
alpha=1;
alpha1=0;
h=0.1;
n=(b-a)/h;
t=linspace(a,b,n+1);
y(1)=alpha;
y1(1)=alpha1;
for i=1:n
    k1=h*f(t(i),y(i),y1(i));
    p1=h*f1(t(i),y(i),y1(i));
    k2=h*f(t(i)+h/2,y(i)+k1/2,y1(i)+p1/2);
    p2=h*f1(t(i)+h/2,y(i)+k1/2,y1(i)+p1/2);
    k3=h*f(t(i)+h/2,y(i)+k2/2,y1(i)+p2/2);
    p3=h*f1(t(i)+h/2,y(i)+k2/2,y1(i)+p2/2);
    k4=h*f(t(i)+h,y(i)+k3,y1(i)+p3);
    p4=h*f1(t(i)+h,y(i)+k3,y1(i)+p3);
    y(i+1)=y(i)+(k1+2*k2+2*k3+k4)/6;
    y1(i+1)=y1(i)+(p1+2*p2+2*p3+p4)/6;
end
plot(t,y,'-b')
hold on
ey=(7/4)*t+(1/2)*(t.^3).*log(t)-(3/4)*(t.^3);
stem (t,ey,'*r')
hold on
legend('approximate solution','actual solution');
fprintf(' t approximate exact error \n');
for i=1:n+1
    fprintf(' %11.5f %11.5f %11.5f %11.5f \n',t(i),y(i),ey(i),abs(ey(i)-y(i)));
end
function z=f(t,y,y1)
z=y1;
function z1=f1(t,y,y1)
z1=(2/t)*y1-(2/(t.^2))*y+t*log(t);