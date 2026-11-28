    % Approximate the solution of the one dimensional heat equation
    %
    % du/dt = d^2u/dx^2,   0 < x < 1,  t > 0
    %
    % which is subject to the initial condition
    % u(x,0) = sin(pi*x),   0 < x < 1
    %
    % and the boundary conditions
    % u(0,t) = u(1,t) = 0,   t >= 0
    %
    % and compare the result at t = 0.04 to the actual solution
    % u(x,t) = exp(-pi^2*t)*sin(pi*x)

    function heat1
    clc
    clf
    clear all
    c=1;%Diffusivity
    a=0;%Left boundary
    b=1;%Right boundary
    tf=0.04;%Final time
    n=10;
    m=12;
    C1=0;%Left Boundary Condition
    C2=0;%Right Boundary Condition
    h=(b-a)/(n-1);%Grid spacing
    k=(tf-0)/(m-1);%Time stepping
    r=(c^2*k)/(h^2)%Stability criterion 0<r<0.5
    u=zeros(n,m);
    for j=1:m
        u(1,j)=C1;
        u(n,j)=C2;
    end
    for i=2:n-1
        u(i,1)=f(a+(i-1)*h);
    end
    for j=2:m
        for i=2:n-1
            u(i,j)=r*u(i-1,j-1)+(1-2*r)*u(i,j-1)+r*u(i+1,j-1);
        end
    end
    disp('Result=')
    num2str(u')
    figure (1)
    surfc(u')
    colorbar
    grid on
    xlabel('x-axis');
    ylabel('t-axis');
    zlabel('u(x,t)');
    hold off
    x=linspace(a,b,n);
    figure (2)
    sum=0.0;
    for i=1:n
        v1(i)=u(i,m);%approximate solution depends on the value of m
        v2(i)=exp(-pi^2*tf)*sin(pi*x(i)); %Exact solution at t=0.04
        sum=sum+(v1(i)-v2(i))*(v1(i)-v2(i));
    end
    rmse=sqrt(sum/length(v1))
    plot(x,v1,'-r')
    hold on
    stem(x,v2,'--*b')
    grid on
    hold off
    title('heat equation')
    legend('approximate solution','actual solution');
    function z=f(x)
    z=sin(pi*x);