% Use the non-linear shooting method to approximate the solution of
% the boundary value problem and compare the result with the actual
% solution.
%
% y'' = 2*y^3,   -1 <= x <= 0,
% y(-1) = 1/2,  y(0) = 1/3,  with h = 0.1
%
% Actual solution:
% y(x) = 1/(x + 3)

function non_linear1
clc
clf
clear all

a = -1;
b = 0;
alpha = 1/2;
beta = 1/3;
h = 0.1;
tol = 0.00001;
m = 10;
k = 1;
n = (b - a)/h;
tk = (beta - alpha)/(b - a);
x = linspace(a, b, n+1);

while k <= m
    u(1) = alpha;
    u1(1) = tk;
    v(1) = 0;
    v1(1) = 1;
    for i = 1:n
        k1 = h*f(x(i), u(i), u1(i));
        k2 = h*f(x(i)+h/2, u(i)+u1(i)*h/2+h*k1/8, u1(i)+k1/2);
        k3 = h*f(x(i)+h/2, u(i)+u1(i)*h/2+h*k1/8, u1(i)+k2/2);
        k4 = h*f(x(i)+h, u(i)+u1(i)*h+h*k3/2, u1(i)+k3);
        u(i+1) = u(i) + h*(u1(i) + (k1+k2+k3)/6);
        u1(i+1) = u1(i) + (k1+2*k2+2*k3+k4)/6;

        p1 = h*g(x(i), v(i), v1(i), u(i), u1(i));
        p2 = h*g(x(i)+h/2, v(i)+v1(i)*h/2+h*p1/8, v1(i)+p1/2, u(i), u1(i));
        p3 = h*g(x(i)+h/2, v(i)+v1(i)*h/2+h*p1/8, v1(i)+p2/2, u(i), u1(i));
        p4 = h*g(x(i)+h, v(i)+v1(i)*h+h*p3/2, v1(i)+p3, u(i), u1(i));
        v(i+1) = v(i) + h*(v1(i) + (p1+p2+p3)/6);
        v1(i+1) = v1(i) + (p1+2*p2+2*p3+p4)/6;
    end
    if (abs(u(n+1) - beta) < tol)
        plot(x, u, '-k')
        hold on
        ey = 1./(x + 3);
        stem(x, ey, '*r')
        hold off
        fprintf('      x       approximate      exact         error\n');
        for i = 1:n+1
            fprintf('%11.5f %11.5f %11.5f %11.5f \n', x(i), u(i), ey(i), abs(ey(i)-u(i)));
        end
        legend('numerical result', 'exact result');
        xlabel('x-axis');
        ylabel('y-axis');
        fprintf('number of iteration = %d\n', k);
        break;
    else
        tk = tk - (u(n+1) - beta)/v(n+1);
        k = k + 1;
    end
end

function u2 = f(x, u, u1)
u2 = 2*u^3;

function v2 = g(x, v, v1, u, u1)
v2 = 6*u^2*v;