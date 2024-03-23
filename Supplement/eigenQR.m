A = [1 2 -1; 1 0 1; 4 -4 5];    % not symmetric
A = A*A'                        % symmetric, so SDT applies
% Since SDT applies, we can use QR iteration
% A = [3 -1 -2; 2 0 -2; 2 -1 -1]-eye(3);
[D,Q] = qrEig(A)

% Check solution using SDT
% A must be recovered from Q*D*Q'
% Diagonals of D are eigenvalues, Q are eigenvectors
Q*D*Q'

function [D,Q] = qrEig(A0)
    n = size(A0,1);
    Ak = A0;        % Note that original A serves as our initial matrix
    Q = eye(n);     % Prepare for accumulating products of Qi
    while (1)
        [Qk,Rk] = GSqr(Ak);
        A = Rk*Qk;
        % This is a good stopping condition since we're looking 
        % at getting the eigenvalues of A (stored in diag(A))
        if norm(diag(A)-diag(Ak)) < 1E-10
            break;
        end
        Ak = A;
        Q = Q*Qk;   % Q will store the eigenvectors as its columns
                    % Cumulative Products(Qi) from i = 1 to k
    end
    D = A;
    Q*A; A0*Q;  % For checking diagonalization theorem
end

function [Q,R] = GSqr(A)
    n = size(A,1);
    Q = A;
    R = zeros(n);
                                % Implied u_1 = v_1
    R(1,1) = norm(Q(:,1));      % r_11 = norm(u_1)
    Q(:,1) = Q(:,1)/R(1,1);     % uhat_1 = u_1/r_11
    for j = 2:n
        % Solve r_ij entries and u_j = v_j - sum(rij * uhat_i), i<j
        for i = 1:j-1
            R(i,j) = Q(:,j)'*Q(:,i);            % r_ij = v_j' * uhat_i, i<j
            Q(:,j) = Q(:,j)-R(i,j)*Q(:,i);      % u_j = u_j - rij
                                                % At i = 1, u_j = v_j
        end
        % Finally, get r_jj and uhat_j
        R(j,j) = norm(Q(:,j));
        Q(:,j) = Q(:,j)/R(j,j);
    end
end