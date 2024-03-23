A=[1 -2 1; 4 -2 1; 1 -2 4];
b=[8 11 17]';
n=size(A,1);

A = [A b] %Augmented Matrix
%Elimination
for j = 1:n
    %Find index of largest magnitude
    [~, k] = max(abs(A(j:end,j)));
    k=(j-1)+k(1);
    %ERO 1
    A([j k], j:end)=A([k j], j:end);
    for i=j+1:n
        if abs(A(i,j))<1E-10
            A(i,j)=0;
            continue;
        end
        %ERO 3
        A(i,j:end)=A(i,j:end)-A(i,j)*A(j,j:end);
    end
end
%REF (A|b)
A
%Backward Sub
x = A(:, end); % x(n) is correct
for i=n-1:-1:1
    x(i)=x(i)-A(i, i+1:end-1)*x(i+1:end);
end
x