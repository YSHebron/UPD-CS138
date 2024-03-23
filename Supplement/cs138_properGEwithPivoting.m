% Assumption: Matrix is square and has a unique solution.
A = [1 -2 1; 4 -2 1; 1 -2 4];
b = [8 11 17]';
n = size(A, 1);

A = [A b] % Augmented Matrix
% Elimination
for j = 1:n
    % Find index of largest magnitude
    [~, k] = max(abs(A(j:end,j)));

    % Adjust k so it reflects row index for entire matrix
    k = (j-1)+k;    % ex. if j = 2 and k = 1, k is actually 2

    % ERO1 switch to get row with pivot into working row
    A([j k], j:end) = A([k j], j:end);

    % ERO2 scale to get a unit pivot
    A(j,j:end) = A(j,j:end)/A(j,j);
    % zero the entries under the unit pivot in column j
    for i = j+1:n
        if abs(A(i,j)) < 1E-10
            A(i,j) = 0;
            continue;
        end
        % ERO3 annihilate
        A(i,j:end) = A(i,j:end)-A(i,j)*A(j,j:end);
    end
end
A %REF (A|b)

%Backward Sub
x = A(:, end); % x(n) is correct
for i = n-1:-1:1
    x(i) = x(i)-A(i,i+1:end-1)*x(i+1:end);
    % ex. x3 = b3
    %     x2 = b2 - a3x3 = b2 - [a3]*[x3]'
    %     x1 = b1 - (a2x2 + a3x3) = b1 - [a2 a3]*[x2 x3]'
end

x   % show solution

