A=[1 -2 0 0 0; -2 1 3 0 0; 0 3 1 -1 0; 0 0 -1 1 4; 0 0 0 4 3]
b=[1 1 1 1 1]';
n=size(A,1);
k1=1; k2=1;  % New parameters

% LU Factorization
% Get QREF of A using only ERO3
% Note that non-unit pivot entries will be allowed
for j=1:n
    if j+k2 == n break; end
    for i=j+1:j+k2
        if abs(A(i,j)) < 1E-10 || A(i,j) == 0
            A(i,j)=0;
            continue;
        end
        if k1 == 0
            break;
        end
        % Store scalar used to annihilate A(i,j) into A(i,j)
        A(i,j)=A(i,j)/A(j,j); % we use A to store both L and U
        A(i,j+1:k1)=A(i,j+1:k1)-A(i,j)*A(j,j+1:k1);
        A
    end
end
A

% Get LU Factorization (L and U stored in A)
U=triu(A)
L=eye(n)+A-triu(A)

[L b]
% Forward Sub Ly = b
y=b;    % y(1) is correct
for i=2:n
    if k2 == 0 break; end
    y(i)=b(i)-A(i,i-k2:i-1)*y(i-k2:i-1);
end
y
[U y]
% Backward Sub Ux = y
x=y;   % x(n) is not correct
x(n)=x(n)/A(n,n);
for i=n-1:-1:1
    if k1 == 0 break; end
    x(i)=(y(i)-A(i,i+1:i+k1)*x(i+1:i+k1))/A(i,i);
end
U*x

x % show solution
A=[1 2 0 0 0; 1 2 4 0 0; 0 1 2 4 0; 0 0 1 2 4; 0 0 0 1 2]
A*x