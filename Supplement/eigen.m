A = [1 2 -1; 1 0 1; 4 -4 5] %lambda = 1,2,3
% A = [1 1; -2 4] %lambda = 2,3
[evec, eval] = powerMethod(A)
[evec, eval] = inverseMethod(A)
[evec, eval] = shiftedMethod(A)

% Observe how these methods still work for
% non-symmetric REAL square matrices

% Power method solves for the eigenpair corresponding
% to the most dominant eigenvalue of A
function [v_new, lambda] = powerMethod(A)
    [n,p] = size(A);
    v0 = rand(n,1);     % Choose initial v randomly
    v = v0;
    while (1)
        v_new = A*v;            
        lambda = v'*v_new/norm(v)^2;    % Rayleigh: v'A*v/v'v
        if norm(v_new - lambda*v) < 1E-6
            break;
        end
        v = v_new/norm(v_new, 'inf');
    end
end

% Inverse Power Method solves for the eigenpair
% corresponding to the least dominant eigenvalue of A
% Uses LU Factorization and Forward and Backward Sub
% as intermediate steps.
function [v_new, lambda] = inverseMethod(A)
    [n,p] = size(A);
    v0 = rand(n,1);     % Choose initial v randomly
    v = v0;
    LU = LUFact(A);      % Take LU Factorization of A
    while (1)
        % To solve for v_new = inv(A)v,
        % use the LU factors on Av_new = LUv_new = v.
        % Forward Sub Ly = v, Backward Sub Uv_new = y
        v_new = forwardAndBack(LU,v);
        lambda = (v'*v_new)/norm(v)^2;
        if norm(v_new - lambda*v) < 1E-6
            break;
        end
        v = v_new/norm(v_new, 'inf');
    end
    lambda = 1/lambda;
end

function [v_new, lambda] = shiftedMethod(A)
    [n,p] = size(A);
    v0 = rand(n,1);
    v = v0;             % choose initial v randomly
    % choose initial lambda as average of max eval and min eval of A.
    [~,min] = powerMethod(A); [~,max] = inverseMethod(A);
    sigma = (min+max)/2;
    lambda = sigma;
    while (1)
        if norm((A-lambda*eye(n))*v) < 1E-6
            break;
        end
        % Solve for v_new in (A-lambda*I)v_new = v;
        % This is inside loop since lambda changes every iter.
        LU = LUFact(A-lambda*eye(n));
        v_new = forwardAndBack(LU,v);
        v = v_new/norm(v_new, 'inf');
        lambda = v'*A*v/norm(v)^2;
    end
    % no need to shift lambda back because it corrects itself
    % as we approximate it with the correct eigenvector
    % Rearrangement of instructions due to lambda
    % being the updating term
end

function A = LUFact(A)
    [n,p] = size(A);
    % Get Quasi-REF of A using only ERO3
    % Note that non-unit pivot entries allowed in QREF
    for j=1:n
        for i=j+1:n
            if abs(A(i,j)) < 1E-10
                A(i,j)=0;
                continue;
            end
            % Store scalar used to annihilate A(i,j) into A(i,j)
            A(i,j)=A(i,j)/A(j,j); % we use A to store both L and U
            A(i,j+1:end)=A(i,j+1:end)-A(i,j)*A(j,j+1:end);
        end
    end
    % At this point A is LU factorized
end

% Companion for LU factorization
function x = forwardAndBack(LU,b)
    % Forward Sub
    [n,p] = size(LU);
    y = b;  % y(1) is correct
    for i=2:n
        y(i)=b(i)-LU(i,1:i-1)*y(1:i-1);
    end

    % Backward Sub
    x = y;   % x(n) is not correct
    x(n)=x(n)/LU(n,n);
    for i=n-1:-1:1
        x(i)=(y(i)-LU(i,i+1:end)*x(i+1:end))/LU(i,i);
    end

    % Remark: x = y = b (space efficient),
    % but we still write them explicitly
    % for ease of understanding.
end