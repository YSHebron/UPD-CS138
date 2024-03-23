% This is a subroutine of Orthogonal Iteration which requires
% A = QR. To implement orthogonal iteration, we insert this
% into the power method while making it work on V instead of v
% and doing GSqr(V) to get Q, which we use instead of V.

A = [1 2 -1; 1 0 1; 4 -4 5]
[n,p] = size(A);    % A is an nxp matrix
Q = A;              % QR Fact can be applied to rect A matrices
R = zeros(p,p);     % hence we make R like this to ensure
                    % that matrix mults (Q*R) are always legal

% Base cases (initially, vj = uj)
R(1,1) = norm(Q(:,1));  % r1 = norm(u1)
Q(:,1) = Q(:,1)/R(1,1); % u1 = u1/r11

% Note that inside Q: vj → uj → uhatj
% Iterate through columns
for j = 2:p
    % Iterate through rows (above curr. diagonal entry)
    for i = 1:j-1
        % rij = vj'*uhati (Update non-diagonal r's, i<j)
        R(i,j) = Q(:,j)'*Q(:,i);
        % uj = vj - rij*uhati (Convert vj into uj)
        % Subtract all needed projections from vj 
        Q(:,j) = Q(:,j) - R(i,j)*Q(:,i);
    end

    % Prepare for next column
    % rjj = norm(uj)
    R(j,j) = norm(Q(:,j));
    % Normalize uj, produce uhatj
    % Reuse R(j,j) for savings
    Q(:,j) = Q(:,j)/R(j,j);     % uj = uj/rjj
end
Q
R
Q*R