% Use the Crank-Nicolson method to approximate the solution of the heat
% equation
%
% du/dt = d^2u/dx^2,   0 < x < 1,  0 < t,
% u(0,t) = u(1,t) = 0,
% u(x,0) = sin(pi*x)
%
% and compare the result to the actual solution
% u(x,t) = exp(-pi^2*t)*sin(pi*x)

function crank1
clc
clf
clear all
c=1;
a=0;
b=1;
tf=0.04;
n=8;
m=7;
C1=0;
C2=0;
h=(b-a)/(n-1);
k=(tf-0)/(m-1);
r=(c^2*k)/(h^2);
s1=2+(2/r);
s2=(2/r)-2;
u=zeros(n,m);
for j=1:m
    u(1,j)=C1;
    u(n,j)=C2;
end
for i=2:n-1
    u(i,1)=f(a+(i-1)*h);
end
A=zeros(n,n);
A(1,1)=1;
A(n,n)=1;
for i=2:n-1
    A(i,i-1)=-1;
    A(i,i)=s1;
    A(i,i+1)=-1;
end
for j=2:m
    F(1,j-1)=C1;
    F(n,j-1)=C2;
    for i=2:n-1
        F(i,j-1)=u(i-1,j-1)+s2*u(i,j-1)+u(i+1,j-1);
    end
    u(:,j)=inv(A)*F(:,j-1);
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
x=linspace(a, b, n);
plot(x,u)
hold on
grid on
title('heat distribution');
for j=1:m
    plot(x,g(x,(j-1)*k),'*k')
end
function z=f(x)
z=sin(x*pi);
function z=g(x,t)
z=exp(-pi^2*t)*sin(pi*x);