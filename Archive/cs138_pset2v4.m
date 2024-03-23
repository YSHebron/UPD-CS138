% item2.m [Hebron, Yenzy]

% imread produces mxnx3 array (3 for RGB)
% im2gray removes third dimension
I = im2gray(imread('peppers.png')); % 'peppers.png' is a built-in
                                    % sample file, no need to submit
I = im2double(I);
% Compute singular value decomposition of I = U*S*V'
% sigk is the number of significant eigenvalues
[U,S,V,sigk] = svdQR(I);
U*S*V'

% Verify correctness of computed SVD
if norm(I - U*S*V', 'inf') < 1E-6
    fprintf("Success!\n");
else
    fprintf("Failed!!!\n")
end

% Generate comparisons of original image I and compressed image I_ap
% using increasing values of k until sigk. We expect the I_ap using k
% values closer to sigk to be of higher quality.
% For peppers.png, sigk is expected to be 384, we use this as our
% test case and we create the subplots around this expectation.
tcs = [1, 40, 80, 120, 160, 200, 240, 280, 320];
for i = 1:9
    k = tcs(i);
    I_ap = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    title(sprintf('I_ap with k = %d', k));
    subplot(3,3,i); imshow(I_ap);
end

%% SVD using QR Iteration (and QR Factorization)
% Every mxn matrix A factors into A = USV'
% U : orthogonal mxm matrix
% S : diagonal mxn matrix
% V : orthogonal nxn matrix
% Wherein orthogonal implies lin. ind. columns (or eigenvectors).
% U and V can be derived from the following equations:
% AA' = U*D_U*U' and A'A = V*D_V*V', where D_U = D_V = SS' = S'S
% Theorem: AA' and A'A are square, symmetric, and non-defective,
% allowing spectral decomposition via QRIter to be applied to either.
function [U,S,V,sigk] = svdQR(A)
    [m,n] = size(A);
    
    % We apply QR Iteration to either AA' or A'A 
    % which would make them converge to D_U or D_V respectively
    % (essentially the same set of eigenvalues).
    % Note that QR Factorization is a subroutine of QR Iteration, and it
    % computes for the transformation matrices Q and R to be used.
    % Moreover, note that as either AA' or A'A converges to the
    % eigenvalues D, Q also converges to the eigenvectors corr. to
    % those eigenvalues, so we also return it.

    % Here, we choose to solve for D_V and V first,
    % and then we compute for U using the equation AV = US
    % (We already know A, right eigenvectors V, and eigenvalues S
    % by that point).
    [D_V, V] = QRIter(A'*A); % V = right singular vectors of A (nxn)
    
    % Clean up eigenvalues D_V
    % Pick sigk by starting from rightmost column of D_V and decrementing
    % until a significant eigenvalue is encountered.
    % We use sigk immediately to improve performance.
    for sigk = size(D_V,2):-1:1
        if D_V(sigk,sigk) >= 1E-6
            break;
        end
    end
    % S must be mxn
    S = zeros(m,n);
    D_V = sqrt(abs(diag(diag(D_V))));
    for j = 1:sigk
        S(j,j) = D_V(j,j);
    end

    % Use obtained V, S and k to compute U in AV = US
    U = zeros(m);
    for j = 1:sigk
        if S(j,j) < 1E-6
            break;
        end
        U(:,j) = A*V(:,j) / S(j,j);
    end
end

%% QR Iteration
% Computes the eigenvalues of A stored in the diagonal matrix D
% corresponding to the eigenvectors of A eventually stored in Q.
% The motivation here is the Spectral Decomposition Theorem in which
% A = X*D*inv(X) = X*D*X' if A is symmetric (meaning inv(X)=X' because
% the eigenvectors of a symmteric matrix are orthogonal, notes p.18).
function [D,Q] = QRIter(A0)
    n = size(A0,1);
    Ak = A0;        % Note that original A serves as our initial matrix
    Q = eye(n);     % Prepare for accumulating products of Qi
    while (1)
        [Qk,Rk] = QRFact(Ak);
        A = Rk*Qk;  % Equiv. to A = Qk'*Ak*Qk (similarity transform.)
        if norm(diag(A)-diag(Ak)) < 0.1     % use large tol to speed-up
            break;      
        end
        Ak = A;
        Q = Q*Qk;   % Q will store the eigenvectors as its columns
                    % Cumulative Products(Qi) from i = 1 to k
    end
    D = A;          % If A is symmteric, Ak converges to diagonals D of A
end

%% QR Factorization equiv. to GS Orthonormalization
% Produces the QR factorization of A
% A : symmetric and non-defective
% Q : contains orthonormalized columns of A
% R : upper triangular matrix relative to the Gram-Schmidt projections.
function [Q,R] = QRFact(A)
    n = size(A,1);  % A is square, so size(A,2) would also do.
    Q = A;          % Initialize Q as A = [v1, v2, ..., vn].
                    % Eventually, Q = [uhat1, uhat2, ..., uhatn].
    R = zeros(n);   % Initialize R as 0 nxn matrix.
                    % This is what we "factor out" of A
    % Base Case
    % Implied u_1 = v_1
    R(1,1) = norm(Q(:,1));          % r_11 = norm(u_1)
    if R(1,1) >= 1E-6
        Q(:,1) = Q(:,1)/R(1,1);     % uhat_1 = u_1/r_11
    end

    % Iterative Case
    % Note that the column vectors of Q undergoes the following transform:
    % Q: v_j → u_j → uhat_j
    for j = 2:n
        % Solve r_ij entries and u_j = v_j - summ(rij * uhat_i), i<j
        % r_ij are the entries of R above the diagonal
        % u_j, by reusing Q(:,j), accumulates the effects of subtracting
        % R(i,j)*Q(:,i), i<j, repeatedly
        % This emulates summ(rij * uhat_i), i<j.
        for i = 1:j-1
            R(i,j) = Q(:,j)'*Q(:,i);            % r_ij = v_j' * uhat_i, i<j
            Q(:,j) = Q(:,j)-(R(i,j)*Q(:,i));    % u_j = v_j - r_ij*uhat_i
                                                % At i = 1, u_j = v_j
        end
        % Finally, get r_jj and uhat_j
        % Get diagonals of R, giving r_jj
        R(j,j) = norm(Q(:,j));      % r_jj = norm(u_j)
        % Normalize jth vector of Q, giving uhat_j
        Q(:,j) = Q(:,j)/R(j,j);     % uhat_j = u_j / r_jj
    end
end
% Note: This subroutine, which always goes through the entire
% matrix Q and is often called by QRiter, is what slows down the code.

