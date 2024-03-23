% imread produces mxnx3 array (3 for RGB)
% im2gray removes third dimension
I = im2gray(imread('peppers.png'));
I = im2double(I);
I = [1 2; 3 4; 5 6]
[U,S,V] = svdQR(I)
U*S*V'
% k = 2;
% I_ap = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
% subplot(121); imshow(I_ap);
% subplot(122); imshow(I);

% U and V are square orthogonal matrices derived from
% AA' = U*D_U*U' and A'A = V*D_V*V'
% Note that while A is not necessarily square,
% AA' and A'A are square and symmetric, allowing SDT → QRIter
function [U,S,V] = svdQR(A)
    [m,n] = size(A);
    % A'A = V*D_V*V'
    [D_V, V] = qrEig(A'*A); % V = right singular vectors of A (nxn)
    
    % Clean up eigenvalues
    % Pick k by starting from rightmost column of D_V and decrementing
    % until a significant eigenvalue is encountered
    for k = size(D_V,2):-1:1
        if D_V(k,k) >= 1E-6
            break;
        end
        D_V(k,:) = zeros(1,n);        % zero insignificant data
        D_V(:,k) = zeros(n,1);
    end
    for j = 2:k
        D_V(j, 1:j-1) = zeros(1,j-1);
        D_V(1:j-1, j) = zeros(j-1,1);
    end
    
    S = sqrt(abs(D_V(1:m, :)));     % S is (mxn)

    % Use obtained V, S and k to compute U in AV = US
    U = zeros(m);
    for j = 1:k
        U(:,j) = A*V(:,j) / S(j,j);
    end
end

function [D,Q] = qrEig(A0)
    n = size(A0,1);
    Ak = A0;        % Note that original A serves as our initial matrix
    Q = eye(n);     % Prepare for accumulating products of Qi
    while (1)
        [Qk,Rk] = GSqr(Ak);
        A = Rk*Qk;
        if norm(diag(A)-diag(Ak)) < 1E-10
            break;
        end
        Ak = A;
        Q = Q*Qk;   % Q will store the eigenvectors as its columns
                    % Cumulative Products(Qi) from i = 1 to k
    end
    D = A;
end

function [Q,R] = GSqr(A)
    n = size(A,1);
    Q = A;
    R = zeros(n);
                                % Implied u_1 = v_1
    R(1,1) = norm(Q(:,1));      % r_11 = norm(u_1)
    Q(:,1) = Q(:,1)/R(1,1);     % uhat_1 = u_1/r_11
    % Note that the column vectors of Q undergoes the following transform:
    % Q: v_j → u_j → uhat_j
    for j = 2:n
        % Solve r_ij entries and u_j = v_j - sum(rij * uhat_i), i<j
        for i = 1:j-1
            R(i,j) = Q(:,j)'*Q(:,i);            % r_ij = v_j' * uhat_i, i<j
            Q(:,j) = Q(:,j)-(R(i,j)*Q(:,i));      % u_j = u_j - rij
                                                % At i = 1, u_j = v_j
        end
        % Finally, get r_jj and uhat_j
        R(j,j) = norm(Q(:,j));      % r_jj = norm(u_j)
        Q(:,j) = Q(:,j)/R(j,j);     % uhat_j = u_j / r_jj
    end
end

