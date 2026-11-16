    % Approximate the one dimensional wave equation
    %
    % d^2u/dt^2 = d^2u/dx^2,   0 < x < 1,  t > 0
    %
    % which is subject to the initial conditions
    % u(x,0) = sin(pi*x)^3   and   du(x,0)/dt = 0   for all 0 <= x <= 1
    %
    % and the boundary conditions
    % u(0,t) = u(1,t) = 0,   t > 0
    %
    % and compare the result at t = 0.04 to the actual solution
    % u(x,t) = (3/4)*sin(pi*x)*cos(pi*t) - (1/4)*sin(3*pi*x)*cos(3*pi*t)

    function wave3
    clc
    clear all
    clf
    c=1;
    a=0;
    b=1;
    tf=0.04;
    n=9;
    m=10;
    c1=0;
    c2=0;
    h=(b-a)/(n-1);
    k=(tf-0)/(m-1);
    r=(c*k)/h;
    s1=1-r^2;
    s2=2-2*r^2;
    U=zeros(n,m);
    for j=1:m
        U(1,j)=c1;
        U(n,j)=c2;
    end
    for i=2:n-1
        U(i,1)=f(a+(i-1)*h);
        U(i,2)=s1*f(a+h*(i-1))+k*g(a+h*(i-1))+(r^2)/2*(f(a+(i-2)*h)+f(a+i*h));
    end
    for j=3:m
        for i=2:n-1
            U(i,j)=s2*U(i,j-1)-U(i,j-2)+(r^2)*(U(i+1,j-1)+U(i-1,j-1));
        end
    end
    disp('Result=')
    num2str(U')
    figure (1)
    clf
    surfc(U')
    colorbar
    grid on
    xlabel('x-axis');
    ylabel('t-axis');
    zlabel('u(x,t)');
    hold off
    figure(2)
    x(1)=a;
    x(n)=b;
    for i=2:n-1
        x(i)=x(i-1)+h;
    end
    title('wave distribution');
    hold on
    grid on
    xlabel('position(x)');
    ylabel('displacement U');
    sum=0.0;
    for i=1:n
        v1(i)=U(i,m);%approximate solution depends on the value of m
        v2(i)=((3/4)*sin(pi*x(i))*cos(pi*tf))-((1/4)*sin(3*pi*x(i))*cos(3*pi*tf)); %Exact solution at t=0.04
        sum=sum+(v1(i)-v2(i))*(v1(i)-v2(i));
    end
    rmse=sqrt(sum/length(v1))
    plot(x,v1,'-r','LineWidth',1.5)
    hold on
    stem(x,v2,'--','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','blue')
    grid on
    legend('approximate solution','actual solution');
    hold off
    function z=f(x)
    z=((sin(pi*x))^3);
    function b=g(x)
    b=0;