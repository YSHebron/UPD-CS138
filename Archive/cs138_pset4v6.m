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
    N = 12;     % N-2: Number of collocation points
    t = linspace(t1,tn,N)'; % time steps
    
    % Chebyshev collocation points
%     N = 10;
%     t = -cos((0:N-1)*pi/(N-1))';
%     t = t1 + ((tn-t1)/2)*(t+1);

    % Basis Function: Multiquadratic Basis
    phi = @(t,tj) sqrt((t-tj).^2 + sigma^2);
    phi_p = @(t,tj) ((t-tj)+1E-6)./phi(t,tj);  % modified with sigma
    phi_pp = @(t,tj) sigma^2./phi(t,tj).^3;

    % Auxiliary Functions for xi and yi and their derivatives
    xi = @(aj,t,tj) sum(aj.*phi(t,tj));
    xi_p = @(aj,t,tj) sum(aj.*phi_p(t,tj));
    xi_pp = @(aj,t,tj) aj.*phi_pp(t,tj);
    yi = @(bj,t,tj) bj.*phi(t,tj);
    yi_p = @(bj,t,tj) bj.*phi_p(t,tj);
    yi_pp = @(bj,t,tj) bj.*phi_pp(t,tj);

    % Initial Guess (random values would do bcz these are just coeffs)
    % Note, theta(:,1) corr to aj, theta(:,2) corr to bj
    theta_0 = randi([-100 100],N,2);

    F = zeros(N,2);
    J = zeros(N*2,N);

    theta = theta_0;       % set theta to initial guess
    while (1)
        % Construct current Collocation Points F
        % Boundary Nodes
        F(1,1) = theta(:,1).*phi(t(1),t) - alpha_x;
        F(N,:) = theta(:,1).*phi(t(N),t) - beta_x;
        F(N+1,:) = theta(:,2).*phi(t(1),t) - alpha_y;
        F(N*2,:) = theta(:,2).*phi(t(N),t) - beta_y;
        % Collocation Points
        for i = 2:N-1
            F(i,:) = xi_pp(theta(:,1),t(i),t) + ...
                    (c/m).*sqrt(xi_p(theta(:,1),t(i),t).^2 + ...
                    yi_p(theta(:,2),t(i),t).^2).* ...
                    xi_p(theta(:,1),t(i),t);
            F(i+N,:) = yi_pp(theta(:,2),t(i),t) + g + ...
                    (c/m).*sqrt(xi_p(theta(:,1),t(i),t).^2 + ...
                    yi_p(theta(:,2),t(i),t).^2).* ...
                    yi_p(theta(:,2),t(i),t);
        end
        % Initiate row-wise summation of F into F_sum
        temp = sum(F,2);
        F_sum = zeros(N,2);
        F_sum(:,1) = temp(1:N); F_sum(:,2) = temp(N+1:end);
        
        % Construct current Jacobian J
        % Boundary Nodes
        J(1,:) = phi(t(1),t);
        J(N,:) = phi(t(N),t);
        J(N+1,:) = phi(t(1),t);
        J(N*2,:) = phi(t(N),t);
        % Internal Nodes
        for i = 2:N-1
            J(i,:) = phi_pp(t(i),t) + (c/m).* ...
                     (2*xi_p(theta(:,1),t(i),t).^2 + ...
                      yi_p(theta(:,2),t(i),t).^2)./ ...
                      sqrt(xi_p(theta(:,1),t(i),t).^2 + ...
                      yi_p(theta(:,2),t(i),t).^2).* ...
                      phi_p(t(i),t);
            J(i+N,:) = phi_pp(t(i),t) + (c/m).* ...
                      (xi_p(theta(:,1),t(i),t).^2 + ...
                      2*yi_p(theta(:,2),t(i),t).^2)./ ...
                      sqrt(xi_p(theta(:,1),t(i),t).^2 + ...
                      yi_p(theta(:,2),t(i),t).^2).* ...
                      phi_p(t(i),t);
        end
        
        % Solve Delta for x and y separately
        
        Delx = J(1:N, :)\F_sum(:,1);
        Dely = J(N+1:end, :)\F_sum(:,2);
        if norm(Delx, 'inf') < 1E-6 && norm(Dely, 'inf') < 1E-6
            break
        else
            theta(:,1) = theta(:,1) - Delx;
            theta(:,2) = theta(:,2) - Dely;
        end
    end
    
    % Coefficients of xi and yi now stored in theta_j (aj and bj resp)
    a = theta(:,1);
    b = theta(:,2);
    z = (t1:(tn-t1)/50:tn)';
    x = zeros(numel(z),numel(z));
    y = zeros(numel(z),numel(z));
    % Internal Nodes
    for i = 2:N
        x = x+a(j)*phi(z,t(j));
        y = y+b(j)*phi(z,t(j));
    end
    % Boundary Nodes
    x = a(1)*phi(z,t(1));
    y = b(1)*phi(z,t(1));
    plot(x,y);
end