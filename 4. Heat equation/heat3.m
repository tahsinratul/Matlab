    % Approximate the solution of the one dimensional heat equation
    %
    % du/dt = (1/pi^2)*d^2u/dx^2,   0 < x < 1,  t > 0
    %
    % which is subject to the initial condition
    % u(x,0) = cos(pi*(x - 1/2)),   0 <= x <= 1
    %
    % and compare the result at t = 0.4 to the actual solution
    % u(x,t) = exp(-t)*cos(pi*(x - 1/2))

    function heat3
    clc
    clear all
    clf
    c=1/pi;%Diffusivity
    a=0;%Left boundary
    b=1;%Right boundary
    tf=0.4;%Final time
    n=10;
    m=12;
    C1=0;%Left Boundary Condition
    C2=0;%Right Boundary Condition
    h=(b-a)/(n-1);%Grid spacing
    k=(tf-0)/(m-1);%Time stepping
    r=(c*c*k)/(h*h)%Stability criterion 0<r<0.5
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
    clf
    surfc(u')
    colorbar
    grid on
    xlabel('x-axis');
    ylabel('t-axis');
    zlabel('u(x,t)');
    hold off
    x(1)=a; x(n)=b;
    for i=2:n-1
        x(i)=x(i-1)+h;
    end
    figure (2)
    sum=0.0;
    for i=1:n
        v1(i)=u(i,m);%approximate solution depends on the value of m
        v2(i)=exp(-tf)*cos(pi*(x(i)-(1/2))); %Exact solution at t=0.4
        sum=sum+(v1(i)-v2(i))*(v1(i)-v2(i));
    end
    rmse=sqrt(sum/length(v1))
    plot(x,v1,'-r','LineWidth',1.5)
    hold on
    stem(x,v2,'--','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','blue')
    grid on
    hold off
    title('heat equation')
    legend('approximate solution','actual solution');
    function z=f(x)
    z=cos(pi*(x-(1/2)));