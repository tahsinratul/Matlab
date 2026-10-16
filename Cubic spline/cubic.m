    % Find the value of y when x = 2, using natural cubic spline interpolation.
    %
    % x:   0.0     0.6      1.5     1.7     1.9     2.1     2.3      2.6       2.8      3.0        3.6       4.7
    % y: -0.8   -0.34   0.59  0.59  0.23    0.1    0.28    1.03    1.5    1.44    0.74     -0.82

    function cubic
    clc
    clear all
    clf
    x=[0.0 0.6 1.5 1.7 1.9 2.1 2.3 2.6 2.8 3.0 3.6 4.7];
    y=[-0.8 -0.34 0.59 0.59 0.23 0.10 0.28 1.03 1.50 1.44 0.74 -0.82];
    p=[1.6 2 2.5 3.5 4.7 0];
    n=length(x)-1;
    for i=1:n
        k(i)=x(i+1)-x(i);
        m(i)=(y(i+1)-y(i))/k(i);
    end
    for i=1:n-1
        A(i,i)=2*(k(i)+k(i+1));
        b(i)=6*(m(i+1)-m(i));
    end
    for i=2:n-1
        A(i-1,i)=k(i);
        A(i,i-1)=k(i);
    end
    z=(A^(-1))*b' ;
    w(1)=0;
    for i=2:n
        w(i)=z(i-1);
    end
    w(n+1)=0;
    plot(x,y,'ok')
    hold on
    for j=1:length(p)
        for i=1:n
            sx=x(i):.01:x(i+1);
            F=(w(i+1)/(6*k(i)))*(sx-x(i)).^3-(w(i)/(6*k(i)))*(sx-x(i+1)).^3+(y(i+1)/k(i)-w(i+1)*k(i)/6)*(sx-x(i))-(y(i)/k(i)-w(i)*k(i)/6)*(sx-x(i+1));
            plot(sx,F,'r')
            if (p(j)>=x(i)&& p(j)<x(i+1))
                sp(j)=(w(i+1)/(6*k(i)))*(p(j)-x(i))^3-(w(i)/(6*k(i)))*(p(j)-x(i+1))^3+(y(i+1)/k(i)-w(i+1)*k(i)/6)*(p(j)-x(i))-(y(i)/k(i)-w(i)*k(i)/6)*(p(j)-x(i+1));
                stem(p(j),sp(j),'*k');
            elseif(p(j)==x(n+1))
                sp(j)=y(n+1);
                stem(p(j),sp(j),'*k');
            end
        end
    end
    fprintf('interpolating points interpolated values\n');
    for i=1:length(p)
        fprintf( '%11.8f %11.8f \n',p(i),sp(i));
    end
