 function cs138_collocation
    nonlinear_colloc();
end

function nonlinear_colloc()
    % Projectile Motion, coupled SODE
    % Boundary Conditions and Assumptions
    t1 = 0; tn = 10;
    alpha_x = 0; beta_x = 8000;
    alpha_y = 0; beta_y = 0;
    global sigma;      % Try sigma = 1,2
    sigma = 1;
    m = 20; c = 0.00032; g = 9.80665;

    % Uniformly spaced collocation points
    global N;
    N = 10;     % Number of collocation points
    t = linspace(t1,tn,N)'; % time steps
    
    % Chebyshev collocation points
%     N = 10;
%     t = -cos((0:N-1)*pi/(N-1))';
%     t = t1 + ((tn-t1)/2)*(t+1);

    % Basis Function: Multiquadratic Basis
    phi = @(t,tj) sqrt((t-tj)^2 + sigma^2);
    phi_p = @(t,tj) (t-tj)/phi(t,tj);
    phi_pp = @(t,tj) sigma^2/phi(t,tj);

    % Initial Guess
    theta_0 = zeros(N*2,1);
    theta_0(1:N) = (alpha_x*(t-tn)-beta_x*(t-t1))/(t1-tn);
    theta_0(N+1:end) = (alpha_y*(t-tn)-beta_y*(t-t1))/(t1-tn);

    F = zeros(N*2,N);
    J = zeros(N*2,N*2);

    theta = theta_0;
    while (1)
        for j = 1:N
            % Boundary Nodes
            F(1,j) = theta(j)*phi(t(1),t(j)) - alpha_x;
            F(N,j) = theta(j)*phi(t(N),t(j)) - beta_x;
            F(N+1,j) = theta(j+N)*phi(t(1),t(j)) - alpha_y;
            F(N*2,j) = theta(j+N)*phi(t(N),t(j)) - beta_y;
            % Collocation Points
            for i = 2:N-1
                F(i,j) = theta(j)*phi_pp(t(i),t(j)) + (c/m)*sqrt((theta(j)*phi_p(t(i),t(j)))^2 + (theta(j+N)*phi_p(t(i),t(j)))^2) * (theta(j)*phi_p(t(i),t(j)));
                F(i+N,j) = theta(j+N)*phi_pp(t(i),t(j)) + g + (c/m)*sqrt((theta(j)*phi_p(t(i),t(j)))^2 + (theta(j+N)*phi_p(t(i),t(j)))^2) * (theta(j+N)*phi_p(t(i),t(j)));
            end
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