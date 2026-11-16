    % Approximate the solution of the Laplace equation
    %
    % d^2u/dx^2 + d^2u/dy^2 = 0,   0 < x < 1,  0 < y < 1
    %
    % subject to the conditions
    % u(x,0) = 10    and   u(x,1) = 0,   0 <= x <= 1
    % u(0,y) = 0     and   u(1,y) = 0,   0 <= y <= 1
    %
    % using the finite difference algorithm.

    function laplace3
    clc
    clear all
    clf
    n=9;
    m=9;
    U=zeros(n,m);
    for i=1:n
        for j=1:m
            h(i,j)=0.0001;
        end
    end
    for i=1:n
        U(i,1)=10;
        U(i,m)=0;
    end
    for j=1:m-1
        U(1,j)=0;
        U(n,j)=0;
    end
    for k=1:1000
        temp1=U';
        for i=2:n-1
            for j=2:m-1
                U(i,j)=(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1))/4;
            end
        end
        temp2=U';
        if (abs(temp2-temp1))<h
            U'
            k
            break
        end
    end
    figure(1)
    clf
    surfc(U')
    colorbar
    grid on
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('U(x,y)');
    hold off
    title('solution of laplace equation');
    figure(2)
    clf
    x=linspace(0,1,n);
    plot(x,U(:,:),'r')
    grid on
    xlabel('x-axis');
    ylabel('y-axis');
    title('solution of laplace equation');