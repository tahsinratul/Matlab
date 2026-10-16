    % Use the linear shooting algorithm to approximate the solution to the
    % boundary value problem and compare the result to the actual solution.
    %
    % y'' = 2*y' - y + x*exp(x) - x,   0 <= x <= 2,
    % y(0) = 0,  y(2) = -4,  with h = 0.1
    %
    % Actual solution:
    % y(x) = (1/6)*x^3*exp(x) - (5/3)*x*exp(x) + 2*exp(x) - x - 2

    function lsa
    clc
    clf
    clear all
    a=0;
    b=2;
    h=0.1;
    n=(b-a)/h;
    x=linspace(a,b,n+1);
    alpha=0;
    beta=-4;
    u(1)=alpha;
    u1(1)=0;
    v(1)=0;
    v1(1)=1;
    for i=1:n
        k1=h*f(x(i),u(i),u1(i));
        k2=h*f(x(i)+h/2,u(i)+u1(i)*h/2+h*k1/8,u1(i)+k1/2);
        k3=h*f(x(i)+h/2,u(i)+u1(i)*h/2+h*k1/8,u1(i)+k2/2);
        k4=h*f(x(i)+h,u(i)+u1(i)*h+h*k3/2,u1(i)+k3);
        u(i+1)=u(i)+h*(u1(i)+(k1+k2+k3)/6);
        u1(i+1)=u1(i)+(k1+2*k2+2*k3+k4)/6;
        p1=h*g(x(i),v(i),v1(i));
        p2=h*g(x(i)+h/2,v(i)+v1(i)*h/2+h*p1/8,v1(i)+p1/2);
        p3=h*g(x(i)+h/2,v(i)+v1(i)*h/2+h*p1/8,v1(i)+p2/2);
        p4=h*g(x(i)+h,v(i)+v1(i)*h+h*p3/2,v1(i)+p3);
        v(i+1)=v(i)+h*(v1(i)+(p1+p2+p3)/6);
        v1(i+1)=v1(i)+(p1+2*p2+2*p3+p4)/6;
    end
    y=u+((beta-u(n+1))/v(n+1))*v;
    plot(x,y,'-b');
    hold on
    ey=(1/6)*x.^3.*exp(x)-(5/3)*x.*exp(x)+2*exp(x)-x-2;
    stem(x,ey,'*k');
    hold off
    legend('numerical result','exact result')
    xlabel('x axis');
    ylabel('y axis');
    fprintf(' x approximate exact error \n');
    for i=1:n+1
        fprintf(' %11.5f %11.5f %11.5f %11.5f \n',x(i),y(i),ey(i),abs(ey(i)-y(i)));
    end
    function z=f(x,u,u1)
    z=2*u1-u+x*exp(x)-x;
    function z1=g(x,v,v1)
    z1=2*v1-v;
