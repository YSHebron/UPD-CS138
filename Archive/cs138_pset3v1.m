% item3.m [Hebron, Yenzy]
syms n k
alpha = 1.35412; beta = 2.67894;

% eigenfunctions, N1 and N2
y1 = @(x,N) 1+symsum( x.^(3*n)/symprod( (3*k+2)*(3*k+3),k,0,n-1 ),n,0,mod(N,2) );
y2 = @(x,N) x+symsum( x.^(3*n+1)/symprod( (3*k+3)*(3*k+4),k,0,n-1 ),n,0,mod(N,3) );

% Airy Functions
Ai = @(x,N) y1(x,N)./(3^(2/3)*alpha) - y2(x,N)./(3^(1/3)*beta);
Bi = @(x,N) y1(x,N)./(3^(1/6)*alpha) + y2(x,N)./(3^(-1/6)*beta);

% Plotting on real domain
x = -15:0.1:5;      % Note: Use linspace for other intervals.
figure(1);
plot(x, Ai(x,5), '-r', x, Ai(x,10), '-b', x, Ai(x,20), '-g');
title("Ai(x):: red: N=5, blue: N=10, green: N=20");
xlabel("x"); ylabel("Ai(x)");
figure(2);
plot(x, Bi(x,5), '-r', x, Bi(x,10), '-b', x, Bi(x,20), '-g');
title("Bi(x):: red: N=5, blue: N=10, green: N=20");
xlabel("x"); ylabel("Bi(x)");
