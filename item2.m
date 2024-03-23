% item2.m [Hebron, Yenzy]

% imread produces mxnx3 array (3 for RGB)
% im2gray removes third dimension
I = im2gray(imread('peppers.png')); % 'peppers.png' is a built-in
                                    % sample file, no need to submit
I = im2double(I);
% Compute singular value decomposition of I = U*S*V'
% sigk is the number of significant eigenvalues
[U,S,V,sigk] = svdQR(I);
I_ap = U*S*V';

% Verify correctness of computed SVD
if norm(I - I_ap, 'inf') < 1E-6
    fprintf("Success!\n");
else
    fprintf("Failed!!!\n")
    fprintf("May have to change correctness tolerance to match that of QRIter.\n")
    fprintf("QRIter tolerance lowered to speed up saturn.png test case.\n")
end

% Generate comparisons of original image I and compressed image I_ap
% using increasing values of k until sigk. We expect the I_ap using k
% values closer to sigk to be of higher quality.
% For peppers.png, sigk is expected to be 384, we use this as our
% test case and we create the subplots around this expectation.
tcs = [1 5 10 20 40 80 160 320 384];
% tcs = 0:sigk/9:sigk;  % Uncomment for more automated testing.
tcs = round(tcs);
norms = zeros(9,1);
for i = 1:9
    k = tcs(i);
    I_ap = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    figure(1);
    subplot(3,3,i); imshow(I_ap);
    title(sprintf('I_{ap} with k = %d', k));
    norms(i) = norm(I-I_ap, 'inf');
end
figure(2);
plot(tcs, norms, 'o-');
xlabel('k'); ylabel('norm(I - Iap)');

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

    % First, we obtain sigk from the smaller of m and n
    if m <= n
        sigk = m;
    else
        sigk = n;
    end
    
    % We then apply QR Iteration to either AA' or A'A 
    % which would make them converge to D_U or D_V respectively
    % (essentially the same set of eigenvalues).
    % Note that QR Factorization is a subroutine of QR Iteration, and it
    % computes for the transformation matrices Q and R to be used.
    % Moreover, note that as either AA' or A'A converges to the
    % eigenvalues D, Q also converges to the eigenvectors corr. to
    % those eigenvalues, so we also return it from QRIter.

    % Here, we choose to solve for D_V and V first,
    % and then we compute for U using the equation AV = US
    % (We already know A, right eigenvectors V, and eigenvalues S
    % by that point).
    [D_V, V] = QRIter(A'*A);
    % D_V and V are nxn
    
    % Compute S which will be mxn.
    % Simultaneously compute inv(S) for later use.
    % S := principal square roots of D_V
    % S must be mxn
    S = zeros(m,n);
    invS = zeros(m,n);
    temp = sqrt(diag(D_V));
    for i = 1:sigk
        S(i,i) = temp(i);       % sqrt of diagonals of D_V
        invS(i,i) = 1/temp(i);  % reciprocal of those diagonals
    end
    % NOTE: Rectangular matrices, like S, don't have inverses in general.
    % What they do have are right and left inverses. Here we will be
    % using that notion instead, wherein invS' is both the left and 
    % right inverse of S, hence we can treat it as a "proper" inverse here.

    % Use now known A, V, and inv(S) to compute U in U=A*V*inv(S)'
    U = A*V*invS';
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

