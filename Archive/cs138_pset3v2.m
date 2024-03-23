% item3.m [Hebron, Yenzy]
syms n k
alpha = 1.35412; beta = 2.67894;

%% Part 1: Plotting on the real domain.
% We later adjust axis limits to make visible the interesting parts
x = -15:0.1:5;      % Note: Use linspace for other intervals.

% 5 term y (also considering 0 terms of y)
y1 = @(x) 1+symsum( x.^(3*n)/symprod( (3*k+2)*(3*k+3),k,0,n-1 ),n,1,1 );
y2 = @(x) x+symsum( x.^(3*n+1)/symprod( (3*k+3)*(3*k+4),k,0,n-1 ),n,1,1 );
Ai = @(x) y1(x)./(3^(2/3)*alpha) - y2(x)./(3^(1/3)*beta);
Bi = @(x) y1(x)./(3^(1/6)*alpha) + y2(x)./(3^(-1/6)*beta);
figure(1);
subplot(1,3,1);
plot(x,double(Ai(x)),'-r');
axis([-10 5 -10 10]);   
title("Ai(x) with 5 terms of y")
xlabel('x'); ylabel('Ai(x)');
figure(2);
subplot(1,3,1);
plot(x,double(Bi(x)),'-b');
axis([-10 5 -10 10]);   
title("Bi(x) with 5 terms of y")
xlabel('x'); ylabel('Bi(x)');

% 10 term y (also considering 0 terms of y)
y1 = @(x) 1+symsum( x.^(3*n)/symprod( (3*k+2)*(3*k+3),k,0,n-1 ),n,1,3 );
y2 = @(x) x+symsum( x.^(3*n+1)/symprod( (3*k+3)*(3*k+4),k,0,n-1 ),n,1,2 );
Ai = @(x) y1(x)./(3^(2/3)*alpha) - y2(x)./(3^(1/3)*beta);
Bi = @(x) y1(x)./(3^(1/6)*alpha) + y2(x)./(3^(-1/6)*beta);
figure(1);
subplot(1,3,2);
plot(x,double(Ai(x)),'-r');
axis([-10 5 -10 10]);
title("Ai(x) with 10 terms of y")
xlabel('x'); ylabel('Ai(x)');
figure(2);
subplot(1,3,2);
plot(x,double(Bi(x)),'-b');
axis([-10 5 -10 10]);   
title("Bi(x) with 10 terms of y")
xlabel('x'); ylabel('Bi(x)');

% 20 term y (also considering 0 terms of y)
y1 = @(x) 1+symsum( x.^(3*n)/symprod( (3*k+2)*(3*k+3),k,0,n-1 ),n,1,6 );
y2 = @(x) x+symsum( x.^(3*n+1)/symprod( (3*k+3)*(3*k+4),k,0,n-1 ),n,1,6 );
Ai = @(x) y1(x)./(3^(2/3)*alpha) - y2(x)./(3^(1/3)*beta);
Bi = @(x) y1(x)./(3^(1/6)*alpha) + y2(x)./(3^(-1/6)*beta);
figure(1);
subplot(1,3,3);
plot(x,double(Ai(x)),'-r');
axis([-10 5 -10 10]);
title("Ai(x) with 20 terms of y")
xlabel('x'); ylabel('Ai(x)');
figure(2);
subplot(1,3,3);
plot(x,double(Bi(x)),'-b');
axis([-10 5 -10 10]);   
title("Bi(x) with 20 terms of y")
xlabel('x'); ylabel('Bi(x)');

% Supplement: Superimposing the Ai(x) and Bi(x) graphs
% using 50 terms of y (0 terms inclusive)
y1 = @(x) 1+symsum( x.^(3*n)/symprod( (3*k+2)*(3*k+3),k,0,n-1 ),n,1,16 );
y2 = @(x) x+symsum( x.^(3*n+1)/symprod( (3*k+3)*(3*k+4),k,0,n-1 ),n,1,15 );
Ai = @(x) y1(x)./(3^(2/3)*alpha) - y2(x)./(3^(1/3)*beta);
Bi = @(x) y1(x)./(3^(1/6)*alpha) + y2(x)./(3^(-1/6)*beta);
figure(3);
plot(x,double(Ai(x)),'-r'); hold on;
plot(x,double(Bi(x)),'-b'); hold off;
axis([-10 5 -10 10]);
title("Ai(x) and Bi(x) with 50 terms of y")

%% Part 2: Plotting on the complex domain.
% We later adjust axis limits to make visible the interesting parts
step = 0.1;
[A,B] = meshgrid(-5:step:5, -5:step:5);

% 20 term y
y1 = @(a,b) 1+symsum( (a+b.*1i).^(3*n)/symprod( (3*k+2)*(3*k+3),k,0,n-1 ),n,1,6 );
y2 = @(a,b) (a+b.*1i)+symsum( (a+b*1i).^(3*n+1)/symprod( (3*k+3)*(3*k+4),k,0,n-1 ),n,1,6 );
Ai = @(a,b) y1(a,b)./(3^(2/3)*alpha) - y2(a,b)./(3^(1/3)*beta);
Bi = @(a,b) y1(a,b)./(3^(1/6)*alpha) + y2(a,b)./(3^(-1/6)*beta);

% Solve for complex Ai and Bi
Aiz = double(Ai(A,B)); fprintf("Solved Aiz!\n")
Biz = double(Bi(A,B)); fprintf("Solved Biz!\n")

figure(4);
subplot(2,2,1);
contourf(A,B,real(Aiz));
title("Re(Ai(x)) with 20 terms of y")
xlabel('a'); ylabel('b'); zlabel('Re(Ai(a,b))')

subplot(2,2,2);
contourf(A,B,imag(Aiz));
title("Im(Ai(x)) with 20 terms of y")
xlabel('a'); ylabel('b'); zlabel('Re(Ai(a,b))')

subplot(2,2,3);
contourf(A,B,real(Biz));
title("Re(Bi(x)) with 20 terms of y")
xlabel('a'); ylabel('b'); zlabel('Re(Ai(a,b))')

subplot(2,2,4);
contourf(A,B,imag(Biz));
title("Im(Bi(x)) with 20 terms of y")
xlabel('a'); ylabel('b'); zlabel('Re(Ai(a,b))')
