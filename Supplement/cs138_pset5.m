

function I=Newton_Cotes(n, a, b, int)
    if (b-a == 0)
        I=0;
        return;
    end
    h = (b-a)/(n-1);
    x = (a:h:b)'; %test points at h intervals
    % CSR 1/3
    w = h/3*ones(n,1);
    w(2:2:end-1) = 4*w(2:2:end-1);
    w(3:2:end-1) = 2*w(3:2:end-1);
    f = int(x);
    I = dot(w,f);
end