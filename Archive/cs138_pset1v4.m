A = [1 -2 -1 3; 2 -4 1 0; 1 -2 2 -3];
b = [1 5 4]';

% A = [1 -2 2 -3; 2 -4 1 0; 1 -2 -1 3];   % Switch up rows
% b = [4 5 1]';
% 
% % Math 40 Test Cases
% A = [0 0 2; 2 10 -4; 1 6 0];
% b = [4 4 9]';
% 
% A = [1 8 0 0 -4; 0 0 1 0 3; 0 0 0 0 0; 0 0 0 1 5];
% b = [-2 1 0 0]';
% 
% A = [1 1 -3 2; 2 -2 2 0; -1 1 -1 0];
% b = [0 0 0]';
% 
% A = [2 0 -1 0; 0 1 2 0; -1 1 2 -1; 0 -1 0 4];
% b = [0 5 2 7]';
% 
% A = [3 0 6 -15 12; 0 2 0 -4 6; 1 0 3 -8 5];
% b = [0 0 0]';
% 
% A = [2 -4 6; 3 -5 8; 0 1 -4];
% b = [2 4 6]';

% Random test cases
% A = [1 2 -3 4; 1 -2 1 0; 5 -2 -3 8; 0 3 1 2];
% b = [0 0 0 0]';
% 
% A = [-2 1 6; -1 5 3; 4 6 12];
% b = [-1 2 -3]';
% 
% A = [2 3 0 0; 1 1 0 0; 4 6 1 -5; 0 5 3 0];
% b = [5 6 7 8]';
% 
% A = [4 0 0 0; 5 1 2 0; 2 6 0 1; -6 3 1 0];
% b = [3 -2 4 -7]';
% 
% A = [0 0; 0 0; 0 0];
% b = [0 0 0]';
% 
% A = [1 0 0; 0 0 0; 0 0 1];
% b = [1 0 3]';
% 
% A = [0 0 0; 0 1 0; 0 0 0; 0 0 1];
% b = [0 0 0 0]';
% 
% A = [2 1 2; 1 0 1; 2 1 0];
% b = [0 -2 4]';
% 
% A = [1 0; 0 1; 1 1];
% b = [3 2 1]';

% Testing eigenspace for iteratively-solved eigen-problems
% A = [1 2 -1; 1 0 1; 4 -4 5];
% A = A-2*eye(3);
% b = [0 0 0]';

x = gaussjordan(A,b)

% To check, see if residual is within our tolerance 1E-10.
% Note that the residual is also a column vector, and the comparison of
% an array with a scalar is done element-wise, returns a "logical array".
% In Matlab, as long as all entries of a matrix is nonzero,
% it is equivalent to true (satisfies if condition).
% Else, it is false. Hence the comparison abs(b-A*x(:,1)) < 1E-10 works.

% Check base case: free variables == 0
if abs(b - A*x(:,1)) < 1E-10
    fprintf('Base case: xhat is correct.\n')
else
    fprintf('Base case: xhat is incorrect!\n')
end

% Check particular (10000 random test cases)
free = size(x,2) - 1;
check = 1;
for i = 1:10000
    plug = [1 randi([-10000 10000], 1,free)]';
    if abs(b - A*(x*plug)) >= 1E-10
        check = 0;
        break
    end
end

if check == 1
    fprintf('All tests passed!\n')
    fprintf('Particular soln xhat*plug is correct.\n')
else
    fprintf('Failed test case detected!\n')
    fprintf('!! Particular soln xhat*plug is incorrect!\n')
end

function ret = gaussjordan(A,b)
    % Get the RREF of A

    [m,n] = size(A);
    p = 0;         % Current number of leading 1s or pivots.
                   % Becomes rank(A) upon getting REF(A).
    A = [A b];  clear b;
    % Get physical limit k of unit pivots
    if m >= n
        limit = n;
    else
        limit = m;
    end
    
    %% Part 1: Gaussian Elimination with Pivoting
    for j = 1:n
        % If number of leading 1s already at matrix limit, stop.
        if p == limit
            break
        end
        
        % ERO1: Switch row with largest magnitude with active row.
        [~,k] = max(abs(A(p+1:end,j)));
        k = p+k;
        % Active row index is given by pivots+1.
        if p+1 ~= k     % If switching by same row ignore switch
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
        % Active row is now given by just p.
        for i = p+1:m
            if abs(A(i,j)) < 1E-10
                A(i,j) = 0;
                continue
            end
            A(i,j:end) = A(i,j:end)-A(i,j)*A(p,j:end);
        end
    end
    
    clear limit;
    A   % Show REF(A)
    
    %% Part 2: REF(A) to RREF(A) (Gauss-Jordan Reduction)
    % Look for pivot points then annihilate entries above pivot points.
    % Note that due to the pivoting implemented here,
    % the nonzero rows are going to be contiguous on the top part
    % of the matrix.
    pcols = [];     % Columns with unit pivots
    for i = p:-1:1
        for j = 1:1:n
            if A(i,j) == 0
                continue
            end
            % Leading entry (1) found at A(i,j)
            for k = i-1:-1:1
                A(k,j:end) = A(k,j:end) - A(k,j)*A(i,j:end);
            end
            % Record pivot column at pcols
            pcols(size(pcols,2)+1) = j;
            
            break
        end
    end
    
    A   % Show RREF(A)

    % Get fcols AFTER getting RREF(A) for better performance.
    % Note: This algo does not see zero columns as free variables. (*)
    pcols = sort(pcols)
    fcols = [];     % Columns with free variables
    % Start checking at col to the right of a pcol
    for i = 1:p
        % Check if all free variable columns found.
        if size(fcols,2) == n-p
            break
        end
        for j = pcols(i)+1:n
            % Also check: if j is already recorded in fcols, move on
            if A(i,j) == 0 || ismember(j,fcols) == 1
                continue
            end
            % Free variable column found at column j
            fcols(size(fcols,2)+1) = j;
        end
    end
    
    %% Part 3: Extracting ret from RREF(A)
    % Note that p also counts the number of nonzero rows.
    % Sort pcols and fcols in preparation for return extraction.
    fcols = sort(fcols)
    
    % Create placeholder matrix for ret, n-p columns for free variables
    % i.e. basis of nullspace(A), and 1 column for particular solution
    ret = zeros(n,n-p);
    particular = zeros(n,1);

    % Particular Solution (Assign constants to rows corresponding to pcols)
    % This assumes a consistent system s.t. no zero row = nonzero constant
    particular(pcols,1) = A(1:p, end);
    
    % If no free variable or trivial input, skip code for free variables.
    % Trivial input := no pivots found
    % Do NOT use here n-p as # of free vars to avoid bug on zero columns
    % (n-p > 0 when there really are no free variables)
    if size(fcols,2) == 0 || p == 0
        ret = particular;
        return;
    end
    
    % Basis of Nullspace(A)
    for j = 1:n-p
        for i = 1:p
            % Use pcols and fcols to "translate" indices bet ret and A
            ret(pcols(i),j) = -A(i,fcols(j));
        end
        % Insert 1 to row representing free variable
        ret(fcols(j),j) = 1;
    end
    
    % Augment particular solution with basis of nullspace of A
    ret = [particular ret];

    return;
end