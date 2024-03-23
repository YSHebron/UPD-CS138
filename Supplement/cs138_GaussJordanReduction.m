A=[1 -2 1; 4 -2 1; 1 -2 4];
b=[8 11 17]';
n=size(A,1);

A=[A b] % Augmented Matrix
%Elimination Part 1 (Regular Gaussian Elimination)
for j=1:n
    %Pivoting
    [~,k]=max(abs(A(j:end,j)));
    k=(j-1)+k;
    %ERO1
    A([j k], j:end) = A([k j], j:end);
    %ERO2
    A(j,j:end)=A(j,j:end)/A(j,j);
    for i = j+1:n
        if abs(A(i,j))<1E-10
            A(i,j)=0;
            continue
        end
        %ERO3
        A(i,j:end)=A(i,j:end)-A(i,j)*A(j,j:end);
    end
end
A

%Elimination Part 2 (Reduction)
%Start from the last column
for j=n:-1:2
    for i=j-1:-1:1
        if abs(A(i,j))<1E-10
            A(i,j)=0;
            continue;
        end
        % ERO 3
        % Observe how 0 entries do not affect the entries
        % above or below it wrt to EROs, hence we can skip
        % over those entries, as done using the staggered index [j end]
        A(i,[j end])=A(i,[j end])-A(i,j)*A(j,[j end]);
    end
end

A
A(:,end)