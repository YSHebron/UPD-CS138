A=[1 -2 1; 4 -2 1; 1 -2 4]; %A=[4 -2 2; -2 2 -4; 2 -4 11];
b=[8 11 17]'; %b=[4 0 -5]';

n = size(A, 1);
A = [A b];
tol = 0.00001;   % define tolerance

% Get [REF(A)|b*]
for j = 1:n
    % Take the row index k of the entry with the largest
    % absolute value in column j. Adjust k using + (j-1).
    % Note: Use the tilde symbol to represent logical NOT
    % or to suppress specific input or output arguments.
    [~,k] = max(abs(A(j:end,j)));
    %k = k+j-1;
    k = (j-1) + k;

    % Then switch row k with your current row j (ERO1).
    % Row j will now refer to switched row.
    temp = A(j,:);
    A(j,:) = A(k,:);
    A(k,:) = temp;

    % alt: A([j k], j:end) = A([k j], j:end);

    % Turn the pivot entry of row j into 1.
    A(j,j:end) = A(j,j:end)/A(j,j); % ERO2
    % line above is optimized to not deal with 0 entries.
    for l = j+1 : n
        % Insert guard for numerically zero entries
        if abs(A(l,j)) < tol
            A(l,j) = 0;
            continue;
        end

        % Annihilate the entries under
        % the pivot entry of row j
        A(l, j:end) = A(l, j:end) - A(l,j)*A(j,j:end); % ERO3
        
        % The A(k,j) corresponds to the
        % appropriate diagonal entry, and
        % A(k,j)*A(j,:) corresponds to scaling row j.
    end
end

% Now A = [REF(A)|b*]
% Apply backsubstitution
% xn is correct at this point
x = A(:, end)  % Get last row of A (matrix b*)
% Note that for i = start : step : stop, step is optional
% Solving for row i (essentially x_i)
for i = n-1 : -1 : 1
    A(i, i+1:end-1)
    x(i+1:end) 
    % dot product is somehow involved with back sub
    x(i) = A(i, end) - A(i, i+1:end-1)*x(i+1:end);
    % the RHS term - the transposed LHS terms
    % ... starting from after the pivot (1) entry
    % ... that is matrix multiplied with the
    % ... back substituted value
end

x   % show solution