    % Use the Crank-Nicolson method to approximate the solution of the heat
    % equation
    %
    % du/dt = d^2u/dx^2,   0 < x < pi,  0 < t,
    % u(0,t) = u(pi,t) = 0,
    % u(x,0) = sin(x)
    %
    % and compare the result to the actual solution
    % u(x,t) = exp(-t)*sin(x)

    function crank2
    clc
    clear all
    clf
    c=1;
    a=0;
    b=pi;
    tf=0.4;
    n=11;
    m=12;
    c1=0;
    c2=0;
    h=(b-a)/(n-1);
    k=(tf-0)/(m-1);
    r=(c*c*k)/(h*h);
    s1=2+(2/r);
    s2=(2/r)-2;
    U=zeros(n,m);
    for j=1:m
        U(1,j)=c1;
        U(n,j)=c2;
    end
    for i=2:n-1
        U(i,1)=f(a+(i-1)*h);
    end
    A(1,1)=1;
    for j=2:m
        for i=2:n-1
            F(1,j-1)=c1;
            F(i,j-1)=U(i-1,j-1)+s2*U(i,j-1)+U(i+1,j-1);
            F(n,j-1)=c2;
            A(1,1)=1;
            A(i,i-1)=-1;
            A(i,i)=s1;
            A(i,i+1)=-1;
            A(n,n)=1;
        end
        U(:,j)=inv(A)*F(:,j-1);
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
    plot(x,U)
    title('heat distribution');
    hold on
    grid on
    for j=1:m
        plot(x,g(x,(j-1)*k),'*k')
    end
    function z=f(x)
    z=sin(x);
    function z=g(x,t)
    z=exp(-t)*sin(x);