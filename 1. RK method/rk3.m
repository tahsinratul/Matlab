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
    alpha1=1;
    alpha2=0;
    h=0.1;
    n=(b-a)/h;
    t=linspace(a,b,n+1);
    u1(1)=alpha1;
    u2(1)=alpha2;
    for i=1:n
        k1=h*f1(t(i),u1(i),u2(i));
        p1=h*f2(t(i),u1(i),u2(i));
        k2=h*f1(t(i)+h/2,u1(i)+k1/2,u2(i)+p1/2);
        p2=h*f2(t(i)+h/2,u1(i)+k1/2,u2(i)+p1/2);
        k3=h*f1(t(i)+h/2,u1(i)+k2/2,u2(i)+p2/2);
        p3=h*f2(t(i)+h/2,u1(i)+k2/2,u2(i)+p2/2);
        k4=h*f1(t(i)+h,u1(i)+k3,u2(i)+p3);
        p4=h*f2(t(i)+h,u1(i)+k3,u2(i)+p3);
        u1(i+1)=u1(i)+(k1+2*k2+2*k3+k4)/6;
        u2(i+1)=u2(i)+(p1+2*p2+2*p3+p4)/6;
    end
    plot(t,u1,'-b')
    hold on
    eu1=(7/4)*t+(1/2)*(t.^3).*log(t)-(3/4)*(t.^3);
    stem (t,eu1,'*r')
    hold on
    legend('approximate solution','actual solution');
    fprintf(' x approximate exact error \n');
    for i=1:n+1
        fprintf(' %11.5f %11.5f %11.5f %11.5f \n',t(i),u1(i),eu1(i),abs(eu1(i)-u1(i)));
    end
    function z1=f1(t,u1,u2)
    z1=u2;
    function z2=f2(t,u1,u2)
    z2=(2/t)*u2-(2/(t.^2))*u1+t*log(t);

