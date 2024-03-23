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
    N = 10;     % N-2: Number of collocation points
    t = linspace(t1,tn,N)'; % time steps
    
    % Chebyshev collocation points
%     N = 10;
%     t = -cos((0:N-1)*pi/(N-1))';
%     t = t1 + ((tn-t1)/2)*(t+1);

    % Basis Function: Multiquadratic Basis
    phi = @(t,tj) sqrt((t-tj).^2 + sigma^2);
    phi_p = @(t,tj) (t-tj)./phi(t,tj);
    phi_pp = @(t,tj) sigma^2./phi(t,tj);

    % Auxiliary Functions for xi and yi and their derivatives
    xi = @(aj,t,tj) aj.*phi(t,tj);
    xi_p = @(aj,t,tj) aj.*phi_p(t,tj);
    xi_pp = @(aj,t,tj) aj.*phi_pp(t,tj);
    yi = @(bj,t,tj) bj.*phi(t,tj);
    yi_p = @(bj,t,tj) bj.*phi_p(t,tj);
    yi_pp = @(bj,t,tj) bj.*phi_pp(t,tj);

    % Initial Guess (use uniformly spaced values from alpha to beta)
    theta_0 = zeros(N*2,1);
    theta_0(1:N) = (alpha_x*(t-tn)-beta_x*(t-t1))/(t1-tn);
    theta_0(N+1:end) = (alpha_y*(t-tn)-beta_y*(t-t1))/(t1-tn);

    F = zeros(N*2,N);       % to be collapsed into F_sum (N*2 x 1)
    J = zeros(N*2,N);

    theta = theta_0;
    while (1)
        % Construct current Collocation Points F
        % Boundary Nodes
        F(1,:) = theta(1:N).*phi(t(1),t) - alpha_x;
        F(N,:) = theta(1:N).*phi(t(N),t) - beta_x;
        F(N+1,:) = theta(N+1:N*2).*phi(t(1),t) - alpha_y;
        F(N*2,:) = theta(N+1:N*2).*phi(t(N),t) - beta_y;
        % Collocation Points
        for i = 2:N-1
            F(i,:) = theta(1:N).*phi_pp(t(i),t) + (c/m).*sqrt((theta(1:N).*phi_p(t(i),t)).^2 + (theta(N+1:N*2).*phi_p(t(i),t)).^2) .* (theta(1:N).*phi_p(t(i),t));
            F(i+N,:) = theta(N+1:N*2).*phi_pp(t(i),t) + g + (c/m).*sqrt((theta(1:N).*phi_p(t(i),t)).^2 + (theta(N+1:N*2).*phi_p(t(i),t)).^2) .* (theta(N+1:N*2).*phi_p(t(i),t));
        end
        
        % Construct current Jacobian J
        % Boundary Nodes
        J(1,1:N) = phi(t(1),t');
        J(N,N+1:N*2) = phi(t(N),t');
        % Internal Nodes
        for i = 2:N-1
            J(i,1:N) = phi_pp(t(i),t) - (-(c/m)*(2*(theta(1:N).*phi_p(t(i),t)).^2 + (theta(N+1:N*2).*phi_p(t(i),t)).^2)./sqrt((theta(1:N).*phi_p(t(i),t)).^2 + (theta(N+1:N*2).*phi_p(t(i),t)).^2));
            J(i,N+1:N*2) = phi_pp(t(i),t) - (-(c/m)*((theta(1:N).*phi_p(t(i),t)).^2 + 2*(theta(N+1:N*2).*phi_p(t(i),t)).^2)./sqrt((theta(1:N).*phi_p(t(i),t)).^2 + (theta(N+1:N*2).*phi_p(t(i),t)).^2)); 
        end
        
        % Clean up NaNs from J
        for i = 1:N*2
            for j = 1:N
                if isnan(J(i,j))
                    J(i,j) = 0;
                end
            end
        end
        
        Sigma = J\F_sum;
        if norm(Sigma, 'inf') < 1E-6
            break
        else
            theta = theta - Sigma;
        end
    end
end