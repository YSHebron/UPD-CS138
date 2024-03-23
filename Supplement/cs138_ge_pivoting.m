A=[1 -2 1; 4 -2 1; 1 -2 4];
b=[8 11 17]';

n = size(A, 1);
A = [A, b];  

% Get [REF(A)|b*]
for j = 1:n
    % Turn the pivot entry of row j into 1.

    A(j,:) = A(j,:)/A(j,j); % ERO2
    for k = j+1 : n
        % Annihilate the entries under
        % the pivot entry of row j
        A(k, :) = A(k, :) - A(k,j)*A(j,:); % ERO3
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