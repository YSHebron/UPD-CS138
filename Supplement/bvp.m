function bvp()
    %linear_shooting()
    %nonlinear_shooting()
    %linear_fda()
    nonlinear_fda()
end

function linear_shooting()
    %% Example: Linear BVP
    % y" = -(2/(1+x^2))y + (2x/1+x^2))y' + 1
    x0 = 0; xf = 4;
    alpha = 1.25;
    beta = -0.95;
    
    % Linear shooting method for DBC
    % Write IVP in terms of system of FODE
    J = @(x,u)[u(2); (2*x/(1+x^2))*u(2)-(2/(1+x^2))*u(1)+1;...
               u(4); (2*x/(1+x^2))*u(4)-(2/(1+x^2))*u(3)+1];
    N = 20; dx = xf/N;
    x = x0:dx:N*dx;

    % Solve for u numerically on the x points.
    [x,u] = ode45(J,x,[alpha 0 alpha 1]');

    % y = c1y1 + c2y2, beta1 = u(end,1), beta2 = u(end,3)
    % Note that y1'(L) = u(end,2) and y2'(L) = u(end,4)
    % which are not needed.
    y = ((beta-u(end,3))/(u(end,1)-u(end,3)))*u(:,1) + ...
        ((u(end,1)-beta)/(u(end,1)-u(end,3)))*u(:,3);

    plot(x,y,'bo-');
end

function nonlinear_shooting()
    %% Example: Nonlinear BVPP
    % y" = y^3 - yy'

    x0 = 1; xf = 2;
    alpha = 0.5;
    beta = 1/3;

    % Nonlinear shooting method for DBC
    % Write IVP in  terms of SFODE
    J = @(x,u)[u(2); u(1)^3-u(1)*u(2); ...
               u(4); -u(1)*u(4)+(3*u(1)^2-u(2))*u(3)];

    N = 50; dx = (xf-x0)/N;
    x = x0:dx:xf;

    % Analytical Solution
    y_ana = 1./(x+1);
    figure(1); plot(x,y_ana,'b'); hold on;

    z = 2; % initial guess
    while (1)
        [~,u] = ode45(J,x,[alpha z 0 1]');
        delz = (u(end,1)-beta)/u(end,3);

        % Plot u(:,1) := y
        plot(x,u(:,1),'bo'); drawnow;
        if abs(delz) < 1E-5
            break;
        else
            z = z-delz; % Newton step
        end
    end
end

function linear_fda()
    %% Example: Linear BVP
    % y" = -2(/1(1+x^2))y + (2x/(1+x^2))y' + 1
    x0 = 0; xf = 4;
    alpha = 1.25;
    beta = -0.95;

    % Linear FDA method
    N = 20; dx = (xf-x0)/(N-1);
    x = x0 + (0:dx:(N-1)*dx)';
    qx = (2.*x./(1+x.^2));
    px = (-2./(1+x.^2));
    rx = ones(size(x));        % Note to return same size as x

    A = zeros(N,N);     % Becomes a banded tridiag matrix
    % Boundary Nodes
    A(1,1) = 1;
    A(N,N) = 1;
    % Internal Nodes
    A(N+2:N+1:end-N) = -(2+dx^2*px(2:end-1));   % 2nd to 2nd to last diag
    A(2:N+1:end-2*N) = 1+(dx/2)*qx(2:end-1);    % internal nodes
    A(2*N+2:N+1:end) = 1-(dx/2)*qx(2:end-1);
    b = dx^2*rx;
    b(1) = alpha;
    b(end) = beta;

    % [A b]
    y = A\b;    % linear solver
    
    figure(1); plot(x,y,'bo-');
end

function nonlinear_fda()
    %% Example: Nonlinear BVP
    % y" = y^3 - yy'
    x0 = 1; xf = 2;
    alpha = 1.25;
    beta = -0.95;
    
    % Nonlinear FDA method
    N = 20; dx = (xf-x0)/(N-1);
    x = x0+(0:dx:(N-1)*dx)';

    % Initial guess
    y = (alpha*(x-xf)-beta*(x-x0))/(x0-xf);

    F = zeros(N,1);
    J = zeros(N,N);
    while (1)
        % Evaluate function vector at y
        % y" = f(x,y,y') = y^3 - yy'
        % Boundary Nodes 
        F(1) = y(1) - alpha;
        F(end) = y(end) - beta;
        % Internal Nodes
        for i = 2:N-1-1
            F(i) = y(i-1)-2*y(i)+y(i+1)-...
                dx^2*(y(i)^3-y(i)*(y(i+1)-y(i-1))/(2*dx));
        end

        % F
        % Evaluate Jacobian at y
        % Boundary Nodes
        J(1,1) = 1;
        J(N,N) = 1;

        % partial wrt y_{i-1}, partial_f/partial_y' = -y
        J(2:N+1:end-2*N) = 1+(dx/2)*(-y(2:end-1));
        % partial wrt y_{i}, partial_f/partial_y' = -y
        J(N+2:N+1:end-N) = -(2+dx^2*(3*(y(2:end-1).^2)-...
                        (y(3:end)-y(1:end-2))/(2*dx)));
        % partial wrt y_{i+1}, partial_f/partial_y' = -y
        J(2*N+2:N+1:end) = 1-(dx/2)*(-y(2:end-1));

        dely = J\F;     % solve for dely
        if norm(dely, 'inf') < 1E-12
            break;
        else
            y = y-dely;
            plot([x0;x;xf], [alpha;y;beta], 'bo'); drawnow;
        end
    end
end