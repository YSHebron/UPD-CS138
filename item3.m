% item3.m [Hebron, Yenzy]
syms n k
alpha = gamma(2/3); beta = gamma(1/3);

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
a = -5:step:5;
b = -5:step:5;
p = numel(a);
q = numel(b);
X = zeros(p,q);
for j = 1:p
    for k = 1:q
        X(j,k) = a(j) + b(k)*i;
    end
end

% Evaluate eigenfunctions on X
y1 = zeros(p,q);
y2 = zeros(p,q);
for k = 1:numel(X)
    y1(k) = 1;
    y2(k) = X(k);
    for n = 1:16
        temp1 = X(k)^(3*n);
        temp2 = X(k)^(3*n+1);
        for l = 0:n-1
            temp1 = temp1 / ((3*l+2)*(3*l+3));
            temp2 = temp2 / ((3*l+3)*(3*l+4));
        end
        y1(k) = y1(k) + temp1;
        y2(k) = y2(k) + temp2;
    end
end

% 20 term y
Ai = y1/(3^(2/3)*alpha) - y2/(3^(1/3)*beta);
Bi = y1/(3^(1/6)*alpha) + y2/(3^(-1/6)*beta);

figure(4);
subplot(2,2,1);
contour(a,b,real(Ai),-5:0.1:5);
title("Re(Ai(x)) with 51 terms of y")
xlabel('a'); ylabel('b');

subplot(2,2,2);
contour(a,b,imag(Ai),-5:0.1:5);
title("Im(Ai(x)) with 51 terms of y")
xlabel('a'); ylabel('b');

subplot(2,2,3);
contour(a,b,real(Bi),-5:0.1:5);
title("Re(Bi(x)) with 51 terms of y")
xlabel('a'); ylabel('b');

subplot(2,2,4);
contour(a,b,imag(Bi),-5:0.1:5);
title("Im(Bi(x)) with 51 terms of y")
xlabel('a'); ylabel('b');
