n = size(A, 1);
A = [A b];
for j = 1:n
    A(j,:) = A(j,:)/A(j,j);
    for k = j+1 : n
        A(k,:) = A(k,:) - A(k,j)*A(j,:);
    end
end
% A = [REF(A)|b*]
% backsubstitution
x = (A:,end);       % -> x_n is correct
for i = n-1 : -1 : 1
    x(i) = A(i, end) - A(i,i+1:end)*x(i+1:end)
end

% Yung star sa Line 6 and 13 ay unsure

for j = 1:n
    k = argmax{A(l, j)} = {A(j,j)}