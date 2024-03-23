% hardcoded sle
A=[1 -2 1; 4 -2 1; 1 -2 4];
b=[8 11 17]';

n = size(A, 1); %take size of first dimension of A
%recall 
A = [A, b];
