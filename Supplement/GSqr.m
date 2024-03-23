% This is a subroutine of QR Iteration which requires Ak = Rk*Qk.
% To implement, we insert this into the power method while making it
% work on V instead of v and doing GSqr(Ak) to to get Qk and Rk.

A = [1 2 -1; 1 0 1; 4 -4 5] % Find A = QR
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
        Q(:,j) = Q(:,j)-R(i,j)*Q(:,i);      % u_j = u_j - rij
                                            % At i = 1, u_j = v_j
    end
    % Finally, get r_jj and uhat_j
    R(j,j) = norm(Q(:,j));      % r_jj = norm(u_j)
    Q(:,j) = Q(:,j)/R(j,j);     % uhat_j = u_j / r_jj
end
Q
R
Q*R