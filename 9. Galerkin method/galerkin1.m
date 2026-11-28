% Use Galerkin weighted residual method to approximate the solution of the
% initial value problem
%
% dy/dx + y = 1,   0 < x < 1,
% y(0) = 0
%
% and compare with the exact solution
% y(x) = 1 - exp(-x)

function galerkin1
clc
clf
clear all
n = 6;
A = zeros(n,n);
B = zeros(n,1);
for i = 1:n
    B(i) = 1/(i+1);
    for j = 1:n
        A(i,j) = j/(i+j) + 1/(i+j+1);
    end
end
C = A\B;
x = linspace(0,1,20);
yG = zeros(size(x));
for j = 1:n
    yG = yG + C(j).*x.^j;
end
yE = 1 - exp(-x);
figure
plot(x,yE,'*k')
hold on
plot(x,yG,'-r')
xlabel('x')
ylabel('y(x)')
title('Galerkin Weighted Residual Method')
legend('Exact Solution','Galerkin Approximation')
grid on
fprintf(' x approximate exact absolute error \n')
for i=1:length(x)
    fprintf('%11.8f %11.8f %11.8f %11.8f\n', x(i), yG(i), yE(i), abs(yE(i)-yG(i)));
end