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
    phi = @(t,tj) sqrt((t-tj).^2 + sigma^2);
    phi_p = @(t,tj) (t-tj)./phi(t,tj);
    phi_pp = @(t,tj) sigma^2./phi(t,tj);

    % Set of Nonlinear Eqns corr to x
    % Unused in practice, but useful as templates.
    % Rename j as k to avoid Matlab confusion with imag. number
    syms k
    xi = @(aj,ti,tj) symsum(aj*phi(ti,tj), k, 1, N);
    xi_p = @(aj,ti,tj) symsum(aj*phi_p(ti,tj), k, 1, N);
    xi_pp = @(aj,ti,tj) symsum(aj*phi_pp(ti,tj), k, 1, N);
    % Set of Nonlinear Eqns corr to y
    yi = @(bj,ti,tj) symsum(bj*phi(ti,tj), k, 1, N);
    yi_p = @(bj,ti,tj) symsum(bj*phi_p(ti,tj), k, 1, N);
    yi_pp = @(bj,ti,tj) symsum(bj*phi_pp(ti,tj), k, 1, N);

    % Initial Guess
    theta_0 = zeros(N*2,1);
    theta_0(1:N) = (alpha_x*(t-tn)-beta_x*(t-t1))/(t1-tn);
    theta_0(N+1:end) = (alpha_y*(t-tn)-beta_y*(t-t1))/(t1-tn);

    F = zeros(N*2,1);
    J = zeros(N,N*2);

    theta = theta_0;
    while (1)
        % Boundary Nodes corr. to x or aj
        F(1) = symsum( theta(k), k, 1, N) - alpha_x;
        F(N) = symsum(theta(k)*phi(t(N)), k, 1, N) - beta_x;
        % Boundary Nodes corr. to y or bj
        F(N+1) = symsum(theta(k+N)*phi(t(1)), k, 1, N) - alpha_y;
        F(end) = symsum(theta(k+N)*phi(t(N)), k, 1, N) - beta_y;
        % Collocation Points corr to x or aj
        for i = 2:N-1
            F(i) = symsum(theta(k)*phi_pp(t(i),t(k)), k, 1, N) ...
                    + (c/m)*sqrt(symsum(theta(k)*phi_p(ti,tj), k, 1, N)^2 ...
                    + symsum(theta(k+N)*phi_p(ti,tj), k, 1, N)^2) ...
                    *symsum(theta(k)*phi_p(ti,tj), k, 1, N);
        end
        % Collocation Points corr to y or bj
        for i = N+1:N*2
            F(i) = symsum(theta(k+N)*phi_pp(ti,tj), k, 1, N) ...
                    + g + (c/m)*sqrt(symsum(theta(k)*phi_p(ti,tj), k, 1, N)^2 ...
                    + symsum(theta(k+N)*phi_p(ti,tj), k, 1, N)^2) ...
                    * symsum(theta(k+N)*phi_p(ti,tj), k, 1, N);
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