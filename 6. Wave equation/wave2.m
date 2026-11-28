% Approximate the one dimensional wave equation
%
% d^2u/dt^2 = d^2u/dx^2,   0 < x < pi,  t > 0
%
% which is subject to the initial conditions
% u(x,0) = sin(x)   and   du(x,0)/dt = 0   for all 0 <= x <= pi
%
% and the boundary conditions
% u(0,t) = u(pi,t) = 0,   t > 0
%
% and compare the result at t = 0.04 to the actual solution
% u(x,t) = cos(t)*sin(x)

function wave2
clc
clf
clear all
c=1;
a=0;
b=pi;
tf=0.04;
n=7;
m=10;
C1=0;
C2=0;
h=(b-a)/(n-1);
k=(tf-0)/(m-1);
r=(c*k)/h;
s1=1-r^2;
s2=2-2*r^2;
u=zeros(n,m);
for j=1:m
    u(1,j)=C1;
    u(n,j)=C2;
end
for i=2:n-1
    u(i,1)=f(a+(i-1)*h);
    u(i,2)=s1*f(a+h*(i-1))+k*g(a+h*(i-1))+(r^2)/2*(f(a+(i-2)*h)+f(a+i*h));
end
for j=3:m
    for i=2:n-1
        u(i,j)=s2*u(i,j-1)-u(i,j-2)+(r^2)*(u(i+1,j-1)+u(i-1,j-1));
    end
end
disp('Result=')
num2str(u')
figure (1)
surfc(u')
colorbar
grid on
xlabel('x-axis');
ylabel('t-axis');
zlabel('u(x,t)');
hold off
figure(2)
x=linspace(a,b,n);
sum=0.0;
for i=1:n
    v1(i)=u(i,m);%approximate solution depends on the value of m
    v2(i)=cos(tf)*sin(x(i)); %Exact solution at t=0.04
    sum=sum+(v1(i)-v2(i))*(v1(i)-v2(i));
end
rmse=sqrt(sum/n)
plot(x,v1,'-r')
hold on
stem(x,v2,'--*b')
grid on
hold off
title('wave distribution');
legend('approximate solution','actual solution');
function z=f(x)
z=sin(x);
function z=g(x)
z=0;