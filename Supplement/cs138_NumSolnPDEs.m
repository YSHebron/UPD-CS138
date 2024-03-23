function cs138_NumSolnPDEs()
    ftcs();     % Forward in Time, Central in Space
    btcs();     % Backward in Time, Central in Space
    cns();      % Crank-Nicolson Scheme
end

function U = ftcs()
    %% Example: Heat Equation (with constant coefficients) (FTCS)
    % delu/delt = D(del^2u/delx^2), x in [x1,xn], t in [t1, tn], D in real+
    x1 = 0; xn = 10; t1 = 0; tn = 1; D = 2;
    % Dirichlet BCs: u(x1,t)=alpha, u(xn,t)=beta;
    alpha = 5; beta = 0;
    % IC: u(x,0)=f(x)

    % Discretization
    n = 20; % number of points in space and time, total of n*n points
    % Note: Prog only works when nx = nt, issue with plotter when nx != nt
    dx = (xn-x1)/(n-1)            % stepsize in space
    x = x1 + (0:dx:(n-1)*dx);     % steps in space 
    dt = (tn-t1)/(n-1)            % stepsize in time, w/ restriction
    restriction = dx^2/(2*D)
    if dt >= restriction
        fprintf("dt is not within restriction! Terminating.\n")
        return;
    end
    t = t1 + (0:dt:(n-1)*dt);     % steps in time
    
    U = ones(n,n);      % initialize U, the ones are placeholders only
    U(1, :) = alpha;
    U(n, :) = beta;
    U(2:end-1,1) = (alpha-beta)/2; %    set init guess here
    for j = 1:n-1
        for i = 2:n-1
            U(i,j+1) = (D*dt)/dx^2*U(i-1,j) + (1-(2*D*dt)/dx^2)*U(i,j) ...
                         + (D*dt)/dx^2*U(i+1,j);
        end
    end

    U   % approximate temp gradient
    plotter(x,t,U,0);
end

function U = btcs()
    %% Example: Heat Equation (with constant coefficients) (BTCS)
    % delu/delt = D(del^2u/delx^2), x in [x1,xn], t in [t1, tn], D in real+
    x1 = 0; xn = 10; t1 = 0; tn = 1; D = 2;
    % Dirichlet BCs: u(x1,t)=alpha, u(xn,t)=beta;
    alpha = 5; beta = 0;
    % IC: u(x,0)=f(x)

    % Discretization
    n = 20; % number of points in space and time, total of n*n points
    % Note: Prog only works when nx = nt, issue with plotter when nx != nt
    dx = (xn-x1)/(n-1);            % stepsize in space
    x = x1 + (0:dx:(n-1)*dx);      % steps in space
    dt = (tn-t1)/(n-1);            % stepsize in time, no restrictions
    t = t1 + (0:dt:(n-1)*dt);      % steps in time

    % Construct Au = b  % Where b = [alpha, ..., dx^2*u_prev(i), ..., beta]
    A = zeros(n,n);     % Becomes almost a constant tridiagonal matrix
                        % Can use LU fact because A is constant
                        % Zeros are (useful) placeholders
    % Boundary Nodes
    A(1,1) = 1;
    A(n,n) = 1;
    % Internal Nodes
    A(n+2:n+1:end-n) = dx^2+2*D*dt;   % coeff of u_i^j (center)
    A(2:n+1:end-2*n) = -D*dt;         % coeff of u_(i-1)^j (left)
    A(2*n+2:n+1:end) = -D*dt;         % coeff of u_(i+1)^j (right)

    % Construct solution matrix V
    % Will dynamically contain:
    % Knowns u_i^(j-1) or u_prev and
    % Unknowns u_i^j or u (should both be 2D)

    % Initial guess on u_i^1 (use this instead of given IC)
    % Note: Use V as container of both u and u_prev
    %  (use V instead of U to avoid confusion with LU)
    V = ones(n,n);
    % V = dx^2*V;       %% mult dx^2 once column is next to be used (!!)
    V(1,1:end) = alpha;
    V(n,1:end) = beta;
    V(2:end-1,1) = (alpha-beta)/2; % set init guess here
    % I choose init guess as avg. of alpha and beta, no abs. since heat
    % or to be more precise temperature can be negative.

    LU = lufact(A,n);   % Workable since A is square

    % No updating, merely iterating through what we already know.
    % Don't confuse this too much with linear FDA.

    % Solve u^j using u^(j-1) (fun to watch contents of V)
    for j = 2:n
        %u = A\b;           % for verifying result of LU fact: okay

        % Solve LUu = u_prev, let Uu = y.
        % Forward sub Ly = u_prev. Solve for y.
        % Take u_prev from V, then modify to reflect b
        %   i.e. mult internal nodes by dx^2
        u_prev = V(:,j-1);
        u_prev(2:end-1,1) = dx^2*u_prev(2:end-1,1);     %% (!!)
        b = u_prev;         % rename for consistency
        y = b;              % y(1) is correct
        for i=2:n
            y(i)=y(i)-LU(i,1:i-1)*y(1:i-1);
        end
    
        % Backward sub Uu = y. Solve for u.
        u = y;
        u(n) = u(n);        % u(1) is correct (boundary)
        for i=n-1:-1:1
            u(i) = (u(i) - LU(i,i+1:end)*u(i+1:end))/LU(i,i);
        end

        V(:,j) = u;  % insert solution to u^j
    end
    
    U = V   % approximate temp gradient
    plotter(x,t,U,1);
end

function U = cns()
    %% Example: Heat Equation (with constant coefficients) (BTCS)
    % delu/delt = D(del^2u/delx^2), x in [x1,xn], t in [t1, tn], D in real+
    x1 = 0; xn = 10; t1 = 0; tn = 1; D = 2;
    % Dirichlet BCs: u(x1,t)=alpha, u(xn,t)=beta;

    % Discretization
    n = 20; % number of points in space and time, total of n*n points
    % Note: Prog only works when nx = nt, issue with plotter when nx != nt
    dx = (xn-x1)/(n-1);            % stepsize in space
    x = x1 + (0:dx:(n-1)*dx);      % steps in space
    dt = (tn-t1)/(n-1);            % stepsize in time, no restrictions
    t = t1 + (0:dt:(n-1)*dt);      % steps in time

    % U1: from FTCS, U2: from BTCS
    U1 = ftcs();
    U2 = btcs();
    U = (U1 + U2) ./ 2;     % take the average

    U
    plotter(x,t,U,2);
end

function [LU] = lufact(A,n)
    LU = A;
    % Iterate through columns, pivoting at diag entries
    for j=1:n
        % Check for failure point.
        if LU(j,j)==0
            disp('Zero pivot encountered');
            return;
        end
        for i=j+1:n
            % Check for numerically zero entries below LU(j,j)
            if abs(LU(i,j)) < 1E-10
                LU(i,j) = 0;
                continue;
            end
            % Get scalar used to annihilate LU(i,j),
            % store scalar in corresponding LU entry
            LU(i,j) = LU(i,j)/LU(j,j);
            % ERO3 to annihilate L(i,j), row subtraction assumed
            LU(i,j+1:end) = LU(i,j+1:end) - LU(i,j)*LU(j,j+1:end);
        end
    end
end

function plotter(x,y,z,fig)
    % Prep for plotting (Borrowed online)
    % https://www.mathworks.com/matlabcentral/answers/387362-how-do-i-
    % create-a-3-dimensional-surface-from-x-y-z-points
    figure(fig*2+1)
    stem3(x, y, z)
    grid on
    
    xv = linspace(min(x), max(x), 100);
    yv = linspace(min(y), max(y), 100);
    [X,Y] = meshgrid(xv, yv);
    Z = griddata(x,y,z,X,Y);
    
    figure(fig*2+2)
    surf(X, Y, Z);
    grid on
    xlabel('space x'); ylabel('time t'); zlabel('heat u(x,t)')
    set(gca, 'ZLim',[0 100])
    shading interp
    view(2)
end