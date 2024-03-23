A = [1 2 -1; 1 0 1; 4 -4 5];    % eigenvals: 3, 2, 1
[n, p] = size(A);
v0 = rand(n,1); 

% % Power Method (Observe how simple this is)
v=v0    % Not yet normalized
it=0;
while (1)
    it = it+1;
    % Iteratively multiply A to v to get v_new
    v_new = A*v;

    % Rayleigh Quotient
    lambda = v'*v_new/norm(v)^2;    % v'A*v/v'v
    

    % At convergence, we want A*v=lambda*v
    if norm(v_new - lambda*v) < 1E-6
        break;
    end

    % Normalized v_new as update to v, inf-norm is most eff. here
    % Avoids overflow error by scaling down v

    v = v_new/norm(v_new, 'inf');
end
lambda
v
it

% %% Inverse Power Method (IPM)
LU = lufact(A,n);
% We can now perform IPM
v=v0;
it=0;
while (1)
    it = it+1;
    % Solve LUv_new = v
    % Forward sub Ly = v
    y = v;      % y(1) is correct
    for i=2:n
        y(i)=y(i)-LU(i,1:i-1)*y(1:i-1);
    end

    % Backward sub Uv_new=y
    v_new = y;
    v_new(n) = v_new(n)/LU(n,n);    % v_new(n) not yet correct, QREF
    for i=n-1:-1:1
        v_new(i) = (v_new(i) - LU(i,i+1:end)*v_new(i+1:end))/LU(i,i);
    end

    % v_new now solved, update lambda via Rayleigh Quotient
    lambda = v'*v_new/norm(v)^2;
    if norm(v_new - lambda*v) < 1E-6
        break;
    end
    
    % Update v, inf-norm most efficient here
    v = v_new/norm(v_new, 'inf');
end
lambda = 1/lambda
v
it

% %% Shifted Inverse Power Method with Rayleigh Quotient
v = v0;
it = 0;
% In choosing initial guess, you may use Gershgorin or take an
% intermediate value between the dominant and recessive eigenvalues.
% Note that a non-defective matrix A of size n has to have n eigenvals.
lambda = 1.6;   % lambda_0 = sigma
while (1)
    it = it+1;

    % Note that lambda changes every iteration, so LU not very efficient
    % Cond. below is as such because at convergence,
    % lambda=sigma → lambda-sigma = 0
    if norm((A-lambda*eye(n))*v) < 1E-6
        break;
    end

    % Gaussian elimination is more efficient here because
    % it only uses back sub, but for this implementation we
    % reuse LU factorization (FS + BS)
    LU = lufact(A-lambda*eye(n),n);
    
    % Solve LUv_new = v
    % Forward sub Ly = v
    y = v;      % y(1) is correct
    for i=2:n
        y(i)=y(i)-LU(i,1:i-1)*y(1:i-1);
    end

    % Backward sub Uv_new=y
    v_new = y;
    v_new(n) = v_new(n)/LU(n,n);    % v_new(n) not yet correct, QREF
    for i=n-1:-1:1
        v_new(i) = (v_new(i) - LU(i,i+1:end)*v_new(i+1:end))/LU(i,i);
    end
    
    % updating v first before lambda speeds up convergence
    % In gen, preemptively using upd. vals. within iteration
    % speeds up convergence. Only valid after stopping cond. checked.
    v = v_new/norm(v_new, 'inf');
    lambda = v'*A*v/norm(v)^2;      % uses updated v
end
lambda
v
it

% LU Factorization for IPM
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