A=[1 -2 1; 4 -2 1; 1 -2 4];
b=[8 11 17]';
n=size(A,1);

% LU Factorization
% Get QREF of A using only ERO3
% Note that non-unit pivot entries will be allowed
for j=1:n
    for i=j+1:n
        if abs(A(i,j)) < 1E-10
            A(i,j)=0;
            continue;
        end
        % Store scalar used to annihilate A(i,j) into A(i,j)
        A(i,j)=A(i,j)/A(j,j); % we use A to store both L and U
        A(i,j+1:end)=A(i,j+1:end)-A(i,j)*A(j,j+1:end);
    end
end
A

% Get LU Factorization (L and U stored in A)
U=triu(A)
L=eye(n)+A-triu(A)

% Forward Sub Ly = b
y=b;    % y(1) is correct
for i=2:n
    y(i)=b(i)-A(i,1:i-1)*y(1:i-1);
end

% Backward Sub Ux = y
x=y;   % x(n) is not correct
x(n)=x(n)/A(n,n);
for i=n-1:-1:1
    x(i)=(y(i)-A(i,i+1:end)*x(i+1:end))/A(i,i);
end

x % show solution