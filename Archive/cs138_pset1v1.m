A = [1 -2 -1 3; 2 -4 1 0; 1 -2 2 -3];
b = [1 5 4]';

% A = [0 0 2; 2 10 -4; 1 6 0];
% b = [4 4 9]';

% A = [1 1 -3 2; 2 -2 2 0; -1 1 -1 0];
% b = [0 0 0]';

% A = [2 -4 6; 3 -5 8; 0 1 -4];
% b = [2 4 6]';

% A = [1 2 -3 4; 1 -2 1 0; 5 -2 -3 8; 0 3 1 2];
% b = [0 0 0 0]';

% A = [-2 1 6; -1 5 3; 4 6 12];
% b = [-1 2 -3]';

% A = [2 3 0 0; 1 1 0 0; 4 6 1 -5; 0 5 3 0];
% b = [5 6 7 8]';

% A = [4 0 0 0; 5 1 2 0; 2 6 0 1; -6 3 1 0];
% b = [3 -2 4 -7]';

% A = [0 0; 0 0; 0 0];
% b = [0 0 0]';

% A = [1 0 0; 0 0 0; 0 0 1];
% b = [1 0 3]';

% A = [0 0 0; 0 1 0; 0 0 0; 0 0 1];
% b = [0 0 0 0]';

% A = [2 1 2; 1 0 1; 2 1 0];
% b = [0 -2 4]';

% A = [1 0; 0 1; 1 1];
% b = [3 2 1]';

% A = [2 0 -1 0; 0 1 2 0; -1 1 2 -1; 0 -1 0 4];
% b = [0 5 2 7]';

% A = [3 0 6 -15 12; 0 2 0 -4 6; 1 0 3 -8 5];
% b = [0 0 0]';

% Shannen
A = [2 -3 1 7; 2 8 -4 5; 1 3 -3 0; -5 2 3 4];
b = [14; -1; 4; -19];
% ans = [1; -3; -4; 1]

A = [2 4 5 7; 1 2 1 -1; -2 -4 1 11];
b = [-26; -4; -10];
% ans = [2 -2 4; 0 1 0; -6 0 -3; 0 0 1]

A = [2 3 -1 -9; 1 2 1 0; -1 2 3 4];
b = [-16; 0; 8];
% ans = [3 -2; -5 3; 7 -4; 0 1]

A = [-1 5 0 0; -2 5 5 2; -3 -1 3 1; 7 6 5 1];
b = [-8; 9; 3; 30];
% ans = [3; -1; 2; 5]

A = [1 2 -4 -1 0; 1 3 -7 0 -1; 1 0 2 -2 3];
b = [32; 33; 22];
% ans = [6 -2 -5; 9 3 2; 0 1 0; -8 0 -1; 0 0 1]

A = [2 1; -1 -1; 3 4; 3 5];
b = [6; -2; 4; 2];
% ans = [4; -2]


A = [2 1 5; 1 -3 -1; 4 -2 6];
b = [10; -2; 12];
% ans = [4 -2; 2 -1; 0 1]

A = [1 2; -3 -1; -2 1];
b = [-4; -3; -7];
% ans = [2; -3]

A = [1 1; -4 -3; 3 2];
b = [1; -2; 1];
% ans = [-1; 2]

A = [1 0 0; 0 1 0; 0 0 1];
b = [1; 2; 3];
% ans = [1; 2; 3]

A = [1 0 0; 0 1 0; 0 0 1];
b = [1; 2; 0];
% ans = [1; 2; 0]

gaussjordan(A,b)

function ret = gaussjordan(A,b)
    % Get the null space of A

    [m,n] = size(A);
    p = 0;         % Number of pivots.
    A = [A b];
    if m >= n
        limit = n;
    else
        limit = m;
    end

    for j = 1:n
        % If number of unit pivots already at matrix limit, stop.
        if p == limit
            break
        end
        
        % Pivoting
        % ERO1: Switch row with largest magnitude with active row.
        [~,k] = max(abs(A(p+1:end,j)));
        k = p+k;
        % Active row index is given by pivots+1.
        if p+1 ~= k
            A([p+1 k], j:end) = A([k p+1], j:end);
        end

        % ERO2: Normalize active row
        if abs(A(p+1,j))<1E-10
            A(p+1,j)=0;
            continue
        end
        % If a non-numerically zero leading entry in active row is found,
        % then we will have a unit pivot. Increase number of found pivots,
        % then proceed to normalize.
        p = p + 1;
        A(p,j:end) = A(p,j:end)/A(p,j);

        % ERO3: Annihilate entries below leading entry of active row.
        % Active row is now given by pivots+1.
        for i = p+1:m
            if abs(A(i,j)) < 1E-10
                A(i,j) = 0;
                continue
            end
            A(i,j:end) = A(i,j:end)-A(i,j)*A(p,j:end);
        end
    end

    refA = A
    
    % Look for pivot points then annihilate entries above pivot points.
    % Note that due to the pivoting implemented here,
    % the nonzero rows are going to be contiguous on the top part
    % of the matrix.
    pcols = [];     % Use for staggered extraction of ret from A
    for i = p:-1:1
        for j = 1:1:n
            if A(i,j) == 0
                continue;
            end
            % Leading entry (1) found at A(i,j)
            for k = i-1:-1:1
                A(k,j:end) = A(k,j:end) - A(k,j)*A(i,j:end);
            end
            [~,l] = size(pcols);
            pcols(l+1) = j;
            break;
        end
    end
    
    rrefA = A

    % Note that p also counts the number of nonzero rows.
    pcols = sort(pcols);
    
    % Create placeholder matrix for return, n-p columns for free variables
    % i.e. basis of nullspace(A), and +1 column for particular solution
    ret = zeros(n,n-p+1);

    % Particular Solution
    ret(pcols,1) = A(1:p, end);
    
    % If no free variable, skip code for free variables.
    if n-p == 0
        return;
    end
    
    % Basis of Nullspace
    fcols = [];
    l = 0;
    for i = 1:p
        k = 2;
        for j = pcols(i)+1:n
            if A(i,j) == 0
                continue;
            end
            ret(pcols(i), k+l) = -A(i,j);
            k = k + 1;
            [~,h] = size(fcols);
            fcols(h+1) = j;
        end
        l = l+1;
    end
    fcols = sort(unique(fcols));
    for i = 1:n-p
        ret(fcols(i),i+1) = 1;
    end

    return;
end