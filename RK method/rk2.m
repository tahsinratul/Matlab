    % Use the RK method for the system to approximate the solution of the
    % following system of first order differential equations and compare
    % the result to the actual solution.
    %
    % u1' = 3*u1 + 2*u2 - (2*t^2 + 1)*exp(2*t),      0 <= t <= 1,  u1(0) = 1
    % u2' = 4*u1 + u2 + (t^2 + 2*t - 4)*exp(2*t),    0 <= t <= 1,  u2(0) = 1,  h = 0.2
    %
    % Actual solution:
    % u1(t) = (1/3)*exp(5*t) - (1/3)*exp(-t) + exp(2*t)
    % u2(t) = (1/3)*exp(5*t) + (2/3)*exp(-t) + t^2*exp(2*t)

    function rk2
    clc
    clf
    clear all
    a=0;
    b=1;
    alpha1=1;
    alpha2=1;
    h=0.2;
    n=(b-a)/h;
    x=linspace(a,b,n+1);
    u1(1)=alpha1;
    u2(1)=alpha2;
    for i=1:n
        k1=h*f1(x(i),u1(i),u2(i));
        p1=h*f2(x(i),u1(i),u2(i));
        k2=h*f1(x(i)+h/2,u1(i)+k1/2,u2(i)+p1/2);
        p2=h*f2(x(i)+h/2,u1(i)+k1/2,u2(i)+p1/2);
        k3=h*f1(x(i)+h/2,u1(i)+k2/2,u2(i)+p2/2);
        p3=h*f2(x(i)+h/2,u1(i)+k2/2,u2(i)+p2/2);
        k4=h*f1(x(i)+h,u1(i)+k3,u2(i)+p3);
        p4=h*f2(x(i)+h,u1(i)+k3,u2(i)+p3);
        u1(i+1)=u1(i)+(k1+2*k2+2*k3+k4)/6;
        u2(i+1)=u2(i)+(p1+2*p2+2*p3+p4)/6;
    end
    eu1=(1/3)*exp(5*x)-(1/3)*exp(-x)+exp(2*x);
    eu2=(1/3)*exp(5*x)+(2/3)*exp(-x)+(x.^2).*exp(2*x);
    plot(x,u1,'-r')
    hold on
    stem (x,eu1,'*r')
    hold on
    plot(x,u2,'-k')
    hold on
    stem(x,eu2,'*k')
    hold on
    legend('Approximate so;ution','actual solution')
    fprintf(' x u1 u2 eu1 eu2 \n');
    for i=1:n+1
        fprintf(' %11.5f %11.5f %11.5f %11.5f %11.5f\n',x(i),u1(i),u2(i),eu1(i),eu2(i));
    end
    function z1=f1(x,u1,u2)
    z1=3*u1+2*u2-(2*x^2+1)*exp(2*x);
    function z2=f2(x,u1,u2)
    z2=4*u1+u2+(x^2+2*x-4)*exp(2*x);