    % Approximate the solution of the Laplace equation
    %
    % d^2u/dx^2 + d^2u/dy^2 = 0,   1 < x < 2,  0 < y < 1
    %
    % subject to the conditions
    % u(x,0) = 2*ln(x)          and   u(x,1) = ln(x^2 + 1),   1 <= x <= 2
    % u(1,y) = ln(y^2 + 1)      and   u(2,y) = ln(y^2 + 4),   0 <= y <= 1
    %
    % using the finite difference algorithm, and compare the result to the
    % actual solution
    % u(x,y) = ln(x^2 + y^2)

    function laplace1
    clc
    clear all
    clf
    n=10;
    m=10;
    a=1;
    b=2;
    c=0;
    d=1;
    u=zeros(n,m);
    x=linspace(a,b,n);
    y=linspace(c,d,m);
    for i=1:n
        u(i,1)=2*log(x(i));
        u(i,m)=log(x(i)^2+1);
    end
    for j=2:m
        u(1,j)=log(y(j)^2+1);
        u(n,j)=log(y(j)^2+4);
    end
    for i=1:n
        for j=1:m
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
    clf
    surfc(u);
    colorbar
    grid on
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('u(x,y)');
    title('Solution of laplace equation using finite difference method')
    for i=1:n
        for j=1:m
            z(i,j)=f(x(i),y(j));
        end
    end
    figure(2)
    plot(x,z,'*')
    hold on
    plot(x,u(:,:),'r')
    xlabel('X- axis')
    ylabel('Y- axis')
    title('Solution of laplace equation using finite difference method')
    function z=f(x,y)
    z=log(x^2+y^2);