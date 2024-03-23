function cs136()
    %I=Newton_Cotes(100)
    %fprintf("%1.6f\n", I)
    fwdEuler(1E-1);
    RK4(1E-1);
end

function RK4(delt)
    t = 0:delt:3;
    y=t;    % just for initialization
    y(1)=1; % initial value
    for i=1:numel(t)-1
        k1=f(t(i), y(i));       u=y(i)+delt/2*k1;
        k2=f(t(i)+delt/2, u);   u=y(i)+delt/2*k2;
        k3=f(t(i)+delt/2, u);   u=y(i)+delt*k3;
        k4=f(t(i)+delt, u);
        k=(k1+2*k2+2*k3+k4)/6;
        y(i+1)=y(i)+delt*k;
    end
    figure(1);
    plot(t,y,'msq'); hold on;
    y_=t.^2+2*t-exp(t)+2;
    plot(t, y_, 'bla-');
end

function fwdEuler(delt)
    t = 0:delt:3;
    y = t;
    y(1) = 1;
    for i = 1:numel(t)-1
        y(i+1) = y(i)+delt*f(t(i), y(i));
    end
    %Plot solution
    figure(1);
    plot(t, y, 'bo'); hold on;
    %Actual solution
    %y_ = exp(-t);
    %plot(t,y_,'r-');
end

function ret=f(t,y)
    % dy/dt = 5y - 6e^-t; soln y = e^-t
    %ret = 5*y-6*exp(-t);
    % dy/dt = y-t^2; soln y=t.^2+2*t-exp(t)+2
    ret = y-t^2;
end

function I=Newton_Cotes(n)
    a = -20;
    b = 20;     % b = pi/2 for Disc Item 2
    h = (b-a)/(n-1);
    x=(a:h:b)'; % (start, step, stop)
    % ones creates a matrix of ones to initialize weights
    %CTR
    %w = h/2*ones(n,1);
    %w(2:end-1) = 2*w(2:end-1);  %mult by 2 except the start and endpoints
    %CSR
    w = h/3*ones(n,1);
    w(2:2:end-1) = 4*w(2:2:end-1);  %mult by 4 per composite (see form)
    w(3:2:end-1) = 2*w(3:2:end-1);  %mult by 2 per composite (see form)
    f = int(x);
    %f(1) = 1;       %redefine f(0)
    I = w'*f;       % get the dot product of weights and function vals
end

function y=int(x)
    %y = sin(x)./x;                             % Disc Item 1
    %y = 1./sqrt(1-(sind(5))^2*(sin(x)).^2);     % Disc Item 2
    y = (1/sqrt(2*pi*2^2))*exp(-(((x-5).^2)./(2*2^2)));
end

function Gauss_Newton()
    y = table2array(readtable("covid_19_data.csv", 'Range', 'C689:C745'));               %csvread("Hubei.csv");
    x = [1:numel(y)]';
    plot(x,y,'bo'); hold on;
    theta = [67000 10 4]';  %parameter vector (theta_k)
    while (1)
        r_n = error(x,y,theta);     %error vector (r_k)
        J_n = Jacobian(x, theta);   %Jacobian of r (Jr(theta_k))
        theta_n = theta-(J_n'*J_n)\(J_n'*r_n);     %theta_k+1
        % Gauss-Newton iteration
        % norm is for general vector distance
        if norm(theta_n - theta) < 1E-4 % tolerance to 4 decimal places
            fprintf('A=%1.2f\n', theta_n(1));
            fprintf('mu = %1.2f\n', theta_n(2));
            fprintf('sigma = %1.2f\n', theta_n(3));
            break;
        else
            theta = theta_n;
        end
    end
    % plot curve fit using new data points at finer intervals
    x_ = x(1):0.1:x(end);
    y_ = theta_n(1)./(1+exp(-(x_-theta_n(2))/theta_n(3)));
    plot(x_, y_, 'b-');
end

