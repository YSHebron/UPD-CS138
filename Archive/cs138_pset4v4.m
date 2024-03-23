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
    phi_p = @(t,tj) ((t-tj)+sigma^2)./phi(t,tj);  % modified with sigma
    phi_pp = @(t,tj) sigma^2./phi(t,tj).^3;

    % Auxiliary Functions for xi and yi and their derivatives
    xi = @(aj,t,tj) aj.*phi(t,tj);
    xi_p = @(aj,t,tj) aj.*phi_p(t,tj);
    xi_pp = @(aj,t,tj) aj.*phi_pp(t,tj);
    yi = @(bj,t,tj) bj.*phi(t,tj);
    yi_p = @(bj,t,tj) bj.*phi_p(t,tj);
    yi_pp = @(bj,t,tj) bj.*phi_pp(t,tj);

    % Initial Guess (use uniformly spaced values from alpha to beta)
    theta_0 = zeros(N,2);
    theta_0(:,1) = (alpha_x*(t-tn)-beta_x*(t-t1))/(t1-tn);
    theta_0(N+1:end) = (alpha_y*(t-tn)-beta_y*(t-t1))/(t1-tn);
    % Replace 0s from initial guess with 1s to avoid NaNs
    for i = 1:N*2
        if theta_0(i) == 0
            theta_0(i) = 1;
        end
    end
    % Clean up NaNs from J, valid since cause of Nans is

    F = zeros(N*2,N);       % to be collapsed into F_sum (N*2 x 1)
    J = zeros(N*2,N);

    theta = theta_0;       % set theta to initial guess
    while (1)
        % Construct current Collocation Points F
        % Boundary Nodes
        F(1,:) = theta(1:N).*phi(t(1),t) - alpha_x;
        F(N,:) = theta(1:N).*phi(t(N),t) - beta_x;
        F(N+1,:) = theta(N+1:N*2).*phi(t(1),t) - alpha_y;
        F(N*2,:) = theta(N+1:N*2).*phi(t(N),t) - beta_y;
        % Collocation Points
        for i = 2:N-1
            F(i,:) = xi_pp(theta(1:N),t(i),t) + ...
                    (c/m).*sqrt(xi_p(theta(1:N),t(i),t).^2 + ...
                    yi_p(theta(N+1:N*2),t(i),t).^2).* ...
                    xi_p(theta(1:N),t(i),t);
            F(i+N,:) = yi_pp(theta(N+1:N*2),t(i),t) + g + ...
                    (c/m).*sqrt(xi_p(theta(1:N),t(i),t).^2 + ...
                    yi_p(theta(N+1:N*2),t(i),t).^2).* ...
                    yi_p(theta(N+1:N*2),t(i),t);
        end
        % Initiate row-wise summation of F into F_rowsum
        F_sum = sum(F,2);
        
        % Construct current Jacobian J
        % Boundary Nodes
        J(1,:) = phi(t(1),t);
        J(N,:) = phi(t(N),t);
        J(N+1,:) = phi(t(1),t);
        J(N*2,:) = phi(t(N),t);
        % Internal Nodes
        for i = 2:N-1
            J(i,:) = phi_pp(t(i),t) + (c/m).* ...
                     (2*xi_p(theta(1:N),t(i),t).^2 + ...
                      yi_p(theta(N+1:N*2),t(i),t).^2)./ ...
                      sqrt(xi_p(theta(1:N),t(i),t).^2 + ...
                      yi_p(theta(N+1:N*2),t(i),t).^2).* ...
                      phi_p(t(i),t);
            J(i+N,:) = phi_pp(t(i),t) + (c/m).* ...
                      (xi_p(theta(1:N),t(i),t).^2 + ...
                      2*yi_p(theta(N+1:N*2),t(i),t).^2)./ ...
                      sqrt(xi_p(theta(1:N),t(i),t).^2 + ...
                      yi_p(theta(N+1:N*2),t(i),t).^2).* ...
                      phi_p(t(i),t);
        end
        
        % Solve Delta for x and y separately
        
        Delx = J\F_sum;
        Dely = J\F_sum(N+1:N*2);
        if norm(Sigma, 'inf') < 1E-6
            break
        else
            theta = theta - Sigma;
        end
    end
end