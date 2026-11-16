    % Approximate the solution of the Laplace equation
    %
    % d^2u/dx^2 + d^2u/dy^2 = 0,   0 < x < 1,  0 < y < 1
    %
    % subject to the conditions
    % u(x,0) = 0     and   u(x,1) = x,   0 <= x <= 1
    % u(0,y) = 0     and   u(1,y) = y,   0 <= y <= 1
    %
    % using the finite difference algorithm, and compare the result to the
    % actual solution
    % u(x,y) = x*y

    function laplace2
    clc
    clear all
    clf
    n=10;
    m=10;
    a=0;
    b=1;
    c=0;
    d=1;
    u=zeros(n,m);
    x=linspace(a,b,n);
    y=linspace(c,d,m);
    for i=1:m
        u(i,1)=0;
        u(i,n)=x(i);
    end
    for j=2:n
        u(1,j)=0;
        u(m,j)=y(j);
    end
    for i=1:m
        for j=1:n
            t(i,j)=0.0001;
        end
    end
    for k=1:1000
        temp1=u';
        for i=2:n-1
            for j=2:m-1
                u(i,j)=(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1))/4;
            end
        end
        temp2=u';
        if (abs(temp1-temp2)) <t
            u'
            k
            break;
        end
    end
    figure(1)
    surfc(u);
    colorbar
    grid on
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('u(x,y)');
    title('Solution of laplace equation using finite difference method')
    figure(2)
    for i=1:n
        for j=1:m
            z(i,j)=f(x(i),y(j));
        end
    end
    figure(2)
    plot(x,z,'*')
    hold on
    plot(x,u,'r')
    xlabel('X- axis')
    ylabel('Y- axis')
    title('Solution of laplace equation using finite difference method')
    function z1=f(x,y)
    z1=x*y;