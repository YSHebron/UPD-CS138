% item1.m [Hebron, Yenzy]

A = [1 -2 -1 3; 2 -4 1 0; 1 -2 2 -3]
b = [1 5 4]'
% x = [2 0 1 0]'

x = gaussjordan(A,b)    % get solution x

% To check, see if norm of residual is within our tolerance 1E-10.
% Check base case: free variables are set to 0
if norm(b - A*x(:,1), 'inf') < 1E-10
    fprintf('Base case: xhat is correct.\n')
else
    fprintf('Base case: xhat is incorrect!\n')
end

% Check particular: free variables can be nonzero
% Test 10000 random test cases
% This might be very slow for huge matrices, so please comment this out
% if ever the performance impact is very bad.
free = size(x,2) - 1;
check = 1;
for i = 1:10000
    % plug serves as the modifier of x.
    % For Pset example, think of plug as [1 s t]'.
    plug = [ 1 randi([-10000 10000], 1, free) ]';
    if norm(b - A*(x*plug), 'inf') >= 1E-10
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

% Solve for the solution x of Ax=b.
function ret = gaussjordan(A,b)
    % Get the RREF of A

    [m,n] = size(A);
    p = 0;         % Current number of leading 1s or pivots. p = pivots
                   % Becomes rank(A) upon getting REF(A).
                   % rank(A) := # of nonzero rows or leading 1s of REF(A).

    % Augment A with b so all transformations done to A also impacts b.
    A = [A b];  clear b;
    % Get physical limit of unit pivots
    % i.e. max # of nonzero rows or leading 1s (pivot cols) of REF(A).
    if m >= n
        limit = n; % if tall/sqr matrix, limit dictated by # of cols n
    else
        limit = m; % if fat matrix, limit dictated by # of rows m
    end
    
    %% Part 1: A to REF(A) (Gaussian Elimination with Pivoting)
    for j = 1:n
        % If number of leading 1s already at matrix limit, stop.
        if p == limit
            break
        end
        
        % Pivoting, increases stability of the method, and also collects
        % all eventually non-zero rows of REF(A) on top rows of the matrix.
        % ERO1: Switch row with largest magnitude with active row.
        [~,k] = max(abs(A(p+1:end,j)));
        % Adjust k to account for subarray A(p+1:end,j).
        % Bug Fix: We use subarray A(p+1:end,j) to avoid switching with
        % "finished rows" (rows already with pivot column) on the
        % next jth (column) iterations.
        k = p+k;
        % k now contains the index of row with largest magnitude.
        % Active row index is given by p+1 (since p starts at 0).
        if p+1 ~= k     % If switching by same row ignore switch
            A([p+1 k], j:end) = A([k p+1], j:end);
        end
        % Also nice to know that the max() function behaves in such a
        % way that we don't switch when row entries in the subarray have
        % the same magnitude.

        % ERO2: Normalize active row
        % Set to zero numerically zero entry, then continue to next
        % column (i.e. there cannot be a pivot column in this column)
        if abs(A(p+1,j))<1E-10
            A(p+1,j)=0;
            continue
        end
        % If a non-numerically zero leading entry in active row is found,
        % then we will have a unit pivot. Increase number of found pivots,
        % then proceed to normalize.
        p = p + 1;  % Active row is now given by just p.
        A(p,j:end) = A(p,j:end)/A(p,j);

        % ERO3: Annihilate entries below leading entry of active row.
        for i = p+1:m
            if abs(A(i,j)) < 1E-10
                A(i,j) = 0;
                continue
            end
            A(i,j:end) = A(i,j:end)-A(i,j)*A(p,j:end);
        end
    end
    
    clear limit;
    REF = A   % Show REF(A)
    
    %% Part 2: REF(A) to RREF(A) (Gauss-Jordan Reduction)
    % Look for pivot points then annihilate entries above pivot points.
    % Note that due to the pivoting implemented here,
    % the nonzero rows are going to be contiguous on the top part
    % of the matrix.
    pcols = [];     % Columns with unit pivots
    for i = p:-1:1  % Start from highest index nonzero row, step backwards.
        for j = 1:1:n   % Scan left to right to guarantee leading entry.
            if A(i,j) == 0
                continue
            end
            % Leading entry (1) found at A(i,j)
            % Annihilate entries above A(i,j:end), just j:end since
            % zero entries don't have effect on other entries.
            for k = i-1:-1:1
                A(k,j:end) = A(k,j:end) - A(k,j)*A(i,j:end);
            end
            % Record pivot column at pcols
            pcols(size(pcols,2)+1) = j;
            
            break
        end
    end
    
    RREF = A   % Show RREF(A)

    % Get fcols AFTER getting RREF(A) for better performance.
    % Note: This algo does not see zero columns as free variables.
    % This is as expected, since zero columns mean variables set to 0.
    pcols = sort(pcols) % Sort first to also get fcols from left to right
    fcols = [];     % Columns with free variables
    % Start checking at col to the right of a pcol
    for i = 1:p
        % Check if all free variable columns found.
        % n-p := expected number of free variables
        if size(fcols,2) == n-p
            break
        end
        for j = pcols(i)+1:n
            % Also check: if j is already recorded in fcols, move on
            if A(i,j) == 0 || ismember(j,fcols) == 1
                continue
            end
            % Free variable column found at column j, record
            fcols(size(fcols,2)+1) = j;
        end
    end
    % Possible optimization: Just record all j's and
    % filter at the end using unique(). Tried this, got Index Error.
    
    %% Part 3: Extracting ret from RREF(A)
    % Note that p also counts the number of nonzero rows.
    % Sort pcols and fcols in preparation for return extraction.
    fcols = sort(fcols)
    
    % Create placeholder matrix for ret, n-p columns for free variables
    % i.e. basis of nullspace(A), and 1 column for particular solution
    ret = zeros(n,n-p);
    particular = zeros(n,1);

    % Particular Solution
    % Assign constants b to rows corresponding to pcols
    % This assumes a consistent system
    % i.e. no zero row corresponding to a nonzero constant
    particular(pcols,1) = A(1:p, end);
    
    % If no free variable or trivial input, skip code for free variables.
    % Trivial input := no pivots found, i.e. zero matrix
    % Do NOT use below n-p as # of free vars to avoid bug on zero columns
    % (It happens that n-p > 0 when there really are no free variables)
    % Just check if fcols contains anything or if there are pivots.
    if size(fcols,2) == 0 || p == 0
        ret = particular;
        return;
    end
    
    % Basis of Nullspace(A)
    for j = 1:n-p
        for i = 1:p
            % Use pcols and fcols to "translate" indices bet ret and A
            % Get negative values since we are "transposing" them.
            ret(pcols(i),j) = -A(i,fcols(j));
            % Breakdown: A(i,fcols(j)) takes free variable coefficients 
            % from A, with j representing which free variable we are
            % dealing with (e.g. in the Pset example, j=1=s, j=2=t).
            % This value is translated to ret(pcols(i),j) which represents
            % how the transposed free variable j values are captured by the
            % solution variables x, with pcols(i) telling which
            % pivot x_k gets the free variable value.
        end
        % Insert 1 to row representing free variable
        ret(fcols(j),j) = 1;
    end
    
    % Augment particular solution with basis of nullspace of A
    ret = [particular ret];

    return;
end