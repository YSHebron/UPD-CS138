% Hardcoded sle
A = [1 -2 1; 4 -2 1; 1 -2 4]; % x = [1 -2 3]'
b = [8 11 17]';

n = size(A, 1);
A = [A b];  % [A|b]

% Get [REF(A)|b*]
for j = 1:n
    % Turn the pivot entry of row j into 1.
    A(j,:) = A(j,:)/A(j,j); % ERO2
    for i = j+1:n
        % Annihilate entries A(i,j) under
        % the pivot entry of row j
        A(i,:) = A(i,:)-A(i,j)*A(j,:); % ERO3
        % Since A(j,j) = 1, lambda = A(i,j)
        % Else lambda = A(i,j)/A(j,j)
    end
end

% Now A = [REF(A)|b*]
% Back substitution
x = A(:,end);  % Get last row of A (matrix b*), x(n)=b(n) is correct
% Note that for i = start : step : stop, step is optional
for i = n-1:-1:1
    % Recall that A is augmented, hence A(i,i+1:end-1)
    x(i) = A(i,end)-A(i,i+1:end-1)*x(i+1:end);
    % the RHS term - the transposed LHS terms
    % ... starting from after the pivot entry
end

x   % show solution