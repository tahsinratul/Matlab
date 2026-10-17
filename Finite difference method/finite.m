    % Use the finite difference method with h = 0.1 to approximate the
    % solution to the boundary value problem
    %
    % d^2y/dx^2 = 2*dy/dx - y(x) + x*exp(x) - x,   0 <= x <= 2,
    % y(0) = 0,  y(2) = -4
    %
    % and compare the result with the actual solution
    % y(x) = (1/6)*x^3*exp(x) - (5/3)*x*exp(x) + 2*exp(x) - x - 2

    function finite
    clc
    clear all
    clf
    a=0;
    b=2;
    alpha=0;
    beta=-4;
    h=0.1;
    n=(b-a)/h;
    x(1)=a;
    x(n+1)=b;
    y(1)=alpha;
    y(n+1)=beta;
    A=zeros(n+1,n+1);
    A(1,1)=1;
    A(n+1,n+1)=1;
    D(1)=alpha;
    D(n+1)=beta;
    for i=2:n
        x(i)= a+(i-1)*h;
        A(i,i-1)= 1+h/2*p(x(i));
        A(i,i)= -h^2*q(x(i))-2;
        A(i,i+1)= 1-h/2*p(x(i));
        D(i)= h^2*r(x(i));
    end
    y= A\D';
    ey=1/6*x.^3.*exp(x)-5/3*x.*exp(x)+2*exp(x)-x-2;
    diff=ey'-y;
    disp('Finite Difference method')
    disp('Result=')
    fprintf(' x y y-exact diff\n')
    disp(([x',y,ey',diff]))
    plot(x,y)
    hold on
    plot(x,ey,'*k')
    grid on
    legend('Approximate solution','Exact solution');
    title('Finite Difference method')
    function z1=p(x)
    z1=2;
    function z2=q(x)
    z2=-1;
    function z3=r(x)
    z3=x*exp(x)-x;