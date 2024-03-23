function eigen_qr()
    tic();
    % The product of any matrix (square or rectangular)
    % and it's transpose is always symmetric.
    A=[1 2 -1; 1 0 1; 4 -4 5]; A=A*A';
    % Diagonals of S are the eigenvalues, Q are the eigenvectors
    [S, Q] = qrEig(A)
    Q*S*Q'
    toc()
end

% QR Iteration
function [S, Q] = qrEig(A0)
    n=size(A0,1);
    Ak=A0;
    Q=eye(n);
    while (1)
        [Qk, Rk] = GSqr(Ak);
        A=Rk*Qk
        if norm(diag(A)-diag(Ak)) < 1E-10
            break;
        end
        Ak=A;
        Q=Q*Qk
    end
    S=A;
end

% Gram-Schmidt (QR Factorization)
function [Q,R]=GSqr(A)
    [n,p]=size(A)
    Q=A;
    R=zeros(p,p);

    R(1,1)=norm(Q(:,1));
    Q(:,1)=Q(:,1)/R(1,1);
    for j=2:p
        for i=1:j-1
            R(i,j)=Q(:,j)'*Q(:,i);
            Q(:,j)=Q(:,j)-R(i,j)*Q(:,i);
        end
        R(j,j)=norm(Q(:,j));
        Q(:,j)=Q(:,j)/R(j,j);
    end
end