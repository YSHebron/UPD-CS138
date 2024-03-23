function myfirstprogram1()
    test1()
    Gauss_Newton()
end

function test1()
    fprintf("Hello world!\n")
    disp("Hello world!")
end

function Gauss_Newton()
    y = table2array(readtable("covid_19_data.csv", 'Range', 'C689:C745'));               %csvread("Hubei.csv");
    x = [1:numel(y)]';
    x
    y
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

function r=error(x,y,theta) % may also work without "label" r
    r = theta(1) ./ (1+exp(-(x-theta(2))/theta(3))) - y;
end

%operator notes:
% ./ : element-wise right division
% .^ : element-wise power

% defining columns of Jacobian (so uses elem wise ops), solved manually
function J=Jacobian(x,theta)   
    J = [1./(1+exp(-(x-theta(2))/theta(3))) ...
            -theta(1)/theta(3)*exp(-(x-theta(2))/theta(3))./(1+exp(-(x-theta(2))/theta(3))).^2 ...
            -theta(1)/theta(3)^2*(x-theta(2)).*exp(-(x-theta(2))/theta(3))./(1+exp(-(x-theta(2))/theta(3))).^2];
end
    