% Use the classical Rk(4,4) method to approximate the solution of the
% initial value problem dy/dx = -5y + 5t^2 + 2t on 0<=t<=1, y(t) = 1/3 with
% h=0.1 and compare to obtained result with the actual solution y(t) = t^2
% + 1/3e^-5t

function rk1
clc
clf
clear all
a=0;
b=1;
alpha=1/3;
h=0.1;
n=(b-a)/h;
t=linspace(a,b,n+1);
y(1)=alpha;
for i=1:n
    k1=h*f(t(i),y(i));
    k2=h*f(t(i)+h/2,y(i)+k1/2);
    k3=h*f(t(i)+h/2,y(i)+k2/2);
    k4=h*f(t(i)+h,y(i)+k3);
    y(i+1)=y(i)+(k1+2*k2+2*k3+k4)/6;
end
plot(t,y,'-b')
hold on
ey=(t.^2)+(1/3)*(exp(-5*t));
stem(t,ey,'*r')
hold on
fprintf(' t approximate exact abserror\n');
for i=1:n+1
    fprintf('%11.5f %11.5f %11.5f %11.5f\n',t(i),y(i),ey(i),abs(ey(i)-y(i)));
end
function z=f(t,y)
z=-5*y+5*t^2+2*t;