function cs138_collocation
    linear_colloc();
end

function linear_colloc()
    %% Example: Linear BVP
    % y''=(2x/(1+x^2))y'-(2/(1+x^2))y+1
    x1 = 0; xn = 4;
    alpha = 1.25;
    beta = -0.95;

    % Collocation Method for Linear BVP
    % y''=(2x/(1+x^2))y'-(2/(1+x^2))y+1

    % Chebyshev collocation points
    N = 10; x = -cos((0:N-1)*pi/(N-1))';
    x = x1 + ((xn-x1)/2)*(x+1);

    r = @(x) ( ones(size(x)) ); %r(x) = 1
    phi = @(x,xj,sigma) ( sqrt((x-xj).^2+sigma^2) );
    Lphi = @(x,xj,sigma,phi_x)( sigma^2./phi_x.^3 ...
                                 -(2*x./(1+x.^2)).*((x-xj)./phi_x) ...
                                 +(2./(1+x.^2)).*phi_x);

    sigma = 1;
    A = zeros(N,N);
    for j = 1:N
        phi_x=phi(x,x(j),sigma);
        A(1,j) = phi_x(1);
        A(2:end-1,j) = Lphi(x(2:end-1),x(j),sigma,phi_x(2:end-1));
        A(end,j) = phi_x(end);
    end
    b = r(x);
    b(1) = alpha;
    b(end) = beta;

    c = A\b;% coefficients          
    z = (x1:(xn-x1)/50:xn)';
    y = c(1)*phi(z,x(1),sigma);
    for j = 2:N
        y = y+c(j)*phi(z,x(j),sigma);
    end

    y_ana = @(x) 1.25+0.4860896526*x-2.25*x.^2+2*x.*atan(x)-log(1+x.^2)/2+(x.^2)/2+(x.^2).*log(1+x.^2)/2;
    figure(1); plot(z,y,'r+'); hold on;
    plot(z,y_ana(z),'b'); hold off;
end