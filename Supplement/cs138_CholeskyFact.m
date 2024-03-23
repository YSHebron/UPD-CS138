% SLE must be symmetric positive definite
A=[4 -2 2; -2 2 -4; 2 -4 11];
b=[4 0 -5]';
n=size(A,1);

% Cholesky Factorization
L=zeros(n,n);
L(1,1)=sqrt(A(1,1));
% Compute base case: L(i,1)
for i=2:n
    L(i,1)=A(i,1)/L(1,1);
end
for j=2:n
    % Compute diagonal entry L(j,j)
    L(j,j)=A(j,j)-L(j,1:j-1)*L(j,1:j-1)';   % A(j,j)-(sum of squares)
    % Check if sqrt is defined on L(j,j)
    if L(j,j)<=0
        disp('Matrix is not positive definite');
        break;
    else
        L(j,j)=sqrt(L(j,j));
    end

    % Compute non-diagonal entry L(i,j)
    % L(i,1:j-1) : L row entries just before j
    % L(j,1:j-1)' : L column entries just before j
    for i = j+1:n
        L(i,j)=(A(i,j)-L(i,1:j-1)*L(j,1:j-1)')/L(j,j);
    end
end

%Forward Sub
y=b; %y(1) is incorrect
y(1)=y(1)/L(1,1);
for i=2:n
    y(i)=(b(i)-L(i,1:i-1)*y(1:i-1))/L(i,i);
end

% Backward Sub
L=L'
x=y; %x(n) is incorrect
x(n)=x(n)/L(n,n);
for i=n-1:-1:1
    x(i)=(y(i)-L(i,i+1:end)*x(i+1:end))/L(i,i)';
end

x % show solution