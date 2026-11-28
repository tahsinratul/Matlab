% Use the finite difference method with h = 0.5 to approximate the
    % solution to the boundary value problem
    %
    % d^2y/dx^2 + y + 1 = 0,   0 <= x <= 1,
    % y(0) = 0,  y(1) = 0
    %
    % and compare the result with the actual solution
    % y = cos(x) + ((1 - cos(1))/sin(1))*sin(x) - 1

    function finite1
    clc
    clf
    clear all
    a=0;
    b=1;
    alpha=0;
    beta=0;
    h=0.5;
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
    ey=cos(x)+((1-cos(1))/sin(1)).*sin(x)-1;
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
    function z=p(x)
    z=0;
    function z=q(x)
    z=-1;
    function z=r(x)
    z=-1;