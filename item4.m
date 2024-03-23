% item4.m [Hebron, Yenzy]
global sigma;   % Try sigma = 1,2
figure(1);
sigma = 1;
subplot(221);
nonlinear_colloc(1,12);
sigma = 2;
subplot(222);
nonlinear_colloc(1,12);
sigma = 1;
subplot(223);
nonlinear_colloc(2,12);
sigma = 2;
subplot(224);
nonlinear_colloc(2,12);

% pts: which kind of colloc pts will be used
% N-2: number of collocation points 
function nonlinear_colloc(pts, N)
    % Projectile Motion, coupled SODE
    % Boundary Conditions and Assumptions
    t1 = 0; tn = 10;
    alpha_x = 0; beta_x = 8000;
    alpha_y = 0; beta_y = 0;
    global sigma;
    m = 20; c = 0.00032;
    g = 9.80665;
    
    %N = 12; % N-2 Number of collocation points
    if pts == 1
        % Use Uniformly spaced collocation points
        t = linspace(t1,tn,N)'; % time steps
    elseif pts == 2
        % Use Chebyshev collocation points
        t = -cos((0:N-1)*pi/(N-1))';
        t = t1 + ((tn-t1)/2)*(t+1); % time steps
    end

    % Basis Function: Multiquadratic Basis
    phi = @(t,tj) sqrt((t-tj).^2 + sigma^2);
    phi_p = @(t,tj) ((t-tj)+1E-6)./phi(t,tj);  % modified with sigma
    phi_pp = @(t,tj) sigma^2./phi(t,tj).^3;

    % Auxiliary Functions for xi and yi and their derivatives
    xi = @(aj,t,tj) sum(aj.*phi(t,tj));
    xi_p = @(aj,t,tj) sum(aj.*phi_p(t,tj));
    xi_pp = @(aj,t,tj) sum(aj.*phi_pp(t,tj));
    yi = @(bj,t,tj) sum(bj.*phi(t,tj));
    yi_p = @(bj,t,tj) sum(bj.*phi_p(t,tj));
    yi_pp = @(bj,t,tj) sum(bj.*phi_pp(t,tj));

    % Initial Guess (random values would do bcz these are just coeffs)
    % Note, theta(:,1) corr to aj, theta(:,2) corr to bj
    theta_0 = randi([-100 100],N,2);

    F = zeros(N,2);
    J = zeros(N*2,N);

    theta = theta_0;       % set theta to initial guess
    while (1)
        % Construct current Collocation Points F
        % Boundary Nodes
        F(1,1) = xi(theta(:,1),t(1),t) - alpha_x;
        F(N,1) = xi(theta(:,1),t(N),t) - beta_x;
        F(1,2) = yi(theta(:,2),t(1),t) - alpha_y;
        F(N,2) = yi(theta(:,2),t(N),t) - beta_y;
        % Collocation Points
        for i = 2:N-1
            F(i,1) = xi_pp(theta(:,1),t(i),t) + ...
                    (c/m)*sqrt(xi_p(theta(:,1),t(i),t)^2 + ...
                    yi_p(theta(:,2),t(i),t)^2)* ...
                    xi_p(theta(:,1),t(i),t);
            F(i,2) = yi_pp(theta(:,2),t(i),t) + g + ...
                    (c/m)*sqrt(xi_p(theta(:,1),t(i),t)^2 + ...
                    yi_p(theta(:,2),t(i),t)^2)* ...
                    yi_p(theta(:,2),t(i),t);
        end
        
        % Construct current Jacobian J
        % Boundary Nodes
        J(1,:) = phi(t(1),t);
        J(N,:) = phi(t(N),t);
        J(N+1,:) = phi(t(1),t);
        J(N*2,:) = phi(t(N),t);
        % Internal Nodes
        for i = 2:N-1
            J(i,:) = phi_pp(t(i),t) + (c/m)* ...
                     (sqrt(xi_p(theta(:,1),t(i),t)^2 + ...
                      yi_p(theta(:,2),t(i),t)^2) + ...
                      xi_p(theta(:,1),t(i),t)^2 / ...
                      sqrt(xi_p(theta(:,1),t(i),t)^2 + ...
                      yi_p(theta(:,2),t(i),t)^2)).* ...
                      phi_p(t(i),t);
            J(i+N,:) = phi_pp(t(i),t) + (c/m)* ...
                     (sqrt(xi_p(theta(:,1),t(i),t)^2 + ...
                      yi_p(theta(:,2),t(i),t)^2) + ...
                      yi_p(theta(:,2),t(i),t)^2 / ...
                      sqrt(xi_p(theta(:,1),t(i),t)^2 + ...
                      yi_p(theta(:,2),t(i),t)^2)).* ...
                      phi_p(t(i),t);
        end
        
        % Solve Delta for x and y separately
%         Delx = J(1:N,:)\F(:,1);
%         Dely = J(N+1:end,:)\F(:,2);
        % We use GE from prev discussions, valid since Jacobian is square
        Delx = ge(J(1:N,:),F(:,1));
        Dely = ge(J(N+1:end,:),F(:,2));
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

    % Reconstructing the trajectory (numerical on top of analytical)
    z = (t1:(tn-t1)/50:tn)';    % Finer interval than t
    % Analytical solution
    x_ana = 800.*z;
    y_ana = (-1/2)*g.*z.^2 + 5*g.*z;
    plot(x_ana,y_ana,'g'); hold on; % green, superimposed
    % Numerical Solution
    x = a(1)*phi(z,t(1));
    y = b(1)*phi(z,t(1));
    for j = 2:N
        x = x+a(j)*phi(z,t(j));
        y = y+b(j)*phi(z,t(j));
    end
    if pts == 1
        plot(x,y,'r+');
        title(sprintf("%d Uniform Colloc Pts. with sigma = %d", N-2, sigma));
    elseif pts == 2
        plot(x,y,'b+');
        title(sprintf("%d Chebyshev Colloc Pts. with sigma = %d", N-2, sigma));
    end
    xlabel('x'); ylabel('y');
    
    % Limit x and y axis for easy comparison
    axis([0 8000 0 8000]); % use y = 0 to 8000 to make angle visible
    
    % Compute approx initial velocity and direction
    vx_0 = xi_p(theta(:,1),t(1),t)
    vy_0 = yi_p(theta(:,2),t(1),t)
    v_0 = sqrt(vx_0^2 + vy_0^2)         % combine velocity components
    angle = rad2deg(atan(y(2)/x(2)))    % compute angle
end

function x = ge(A,b)
    % Assumption: Matrix is square and has a unique solution.
    n = size(A, 1);
    
    A = [A b]; % Augmented Matrix
    % Elimination
    for j = 1:n
        % Find index of largest magnitude
        [~, k] = max(abs(A(j:end,j)));
    
        % Adjust k so it reflects row index for entire matrix
        k = (j-1)+k;    % ex. if j = 2 and k = 1, k is actually 2
    
        % ERO1 switch to get row with pivot into working row
        A([j k], j:end) = A([k j], j:end);
    
        % ERO2 scale to get a unit pivot
        A(j,j:end) = A(j,j:end)/A(j,j);
        % zero the entries under the unit pivot in column j
        for i = j+1:n
            if abs(A(i,j)) < 1E-10
                A(i,j) = 0;
                continue;
            end
            % ERO3 annihilate
            A(i,j:end) = A(i,j:end)-A(i,j)*A(j,j:end);
        end
    end
    % A %REF (A|b)
    
    %Backward Sub
    x = A(:, end); % x(n) is correct
    for i = n-1:-1:1
        x(i) = x(i)-A(i,i+1:end-1)*x(i+1:end);
        % ex. x3 = b3
        %     x2 = b2 - a3x3 = b2 - [a3]*[x3]'
        %     x1 = b1 - (a2x2 + a3x3) = b1 - [a2 a3]*[x2 x3]'
    end
end


