% Use the Galerkin method to approximate the solution of the boundary
% value problem
%
% d^2y/dx^2 + y + 2*x*(1 - x) = 0,   0 < x < 1,
% y(0) = 0,  y(1) = 1/2
%
% and compare with the exact solution
% y(x) = 4*cos(x) - 2*x + 2*x^2 - ((8*cos(1) - 9)/(2*sin(1)))*sin(x) - 4

function galerkin
clc
clear all
clf
n = 6;
A = zeros(n,n);
B = zeros(n,1);
for i = 1:n
    B(i) = -5/(2*(i+2))+9/(2*(i+3))-2/(i+4);
    for j = 1:n
        A(i,j) = j*(j-1)/(i+j-1)-(2*j^2)/(i+j)+(j^2+j+1)/(i+j+1)-2/(i+j+2)+1/(i+j+3);
    end
end
C = A\B;
x = linspace(0,1,20);
yG = x/2;
for j = 1:n
    yG = yG + C(j).*x.^j.*(1-x);
end
yE = 4*cos(x)-2*x+2*x.^2 -((8*cos(1)-9)/(2*sin(1))).*sin(x)-4;
figure
plot(x,yE,'*k')
hold on
plot(x,yG,'-r')
xlabel('x')
ylabel('y(x)')
title('Galerkin Weighted Residual Method')
legend('Exact Solution','Galerkin Approximation','Location','best')
grid on
fprintf(' x approximate exact absolute error \n')
for i=1:length(x)
    fprintf('%11.8f %11.8f %11.8f %11.8f\n', x(i), yG(i), yE(i), abs(yE(i)-yG(i)));
end