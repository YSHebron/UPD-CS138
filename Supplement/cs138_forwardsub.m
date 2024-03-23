A=[1 0 0; 2 3 0; 4 5 6];
b=[4 5 6];
x=A(:, 1);
n = size(A,1);
for j = 1:n
    if A(j,j) < 1E-10
        break;
    end
    x(j) = b(j)/A(j,j);
    for i = j+1 : n
        b(i) = b(i) - A(i,j) * x(j);   % preemptively updates next b(j)
    end
end
    
x