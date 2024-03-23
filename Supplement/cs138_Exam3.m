function cs138_Exam3()
    heat_eqn();
end

function heat_eqn()
    % Taking 5 terms only
    D = 0.1; e = exp(1);

    u_ana = @(x,t) (- (4/pi^2)*e^(-pi^2*D*t)*cos(pi*x) ...
                    + (1/pi^2)*e^(-4*pi^2*D*t)*cos(2*pi*x) ...
                    - (4/(9*pi^2))*e^(-9*pi^2*D*t)*cos(3*pi*x) ...
                    + (1/(4*pi^2))*e^(-16*pi^2*D*t)*cos(4*pi*x) ...
                    - (4/(25*pi^2))*e^(-25*pi^2*D*t)*cos(5*pi*x));

    fsurf(u_ana, [-1 1 0 1]);
    title("Simple Heat Equation (PBC)")
    xlabel('space x') 
    ylabel('time t')
    zlabel('heat u(x,t)')
    view(2)
end