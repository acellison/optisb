function [uk,u] = rk4(f,h,t0,u0,n)
    if nargout>1
        u = zeros(n,length(u0),class(u0));
    end
    tk = t0;
    uk = u0;
    for i=1:n
        k1 = f(tk,     uk);
        k2 = f(tk+.5*h,uk+.5*h*k1);
        k3 = f(tk+.5*h,uk+.5*h*k2);
        k4 = f(tk+   h,uk+   h*k3);
        uk = uk+(h/6)*(k1+2*k2+2*k3+k4);
        tk = tk+h;
        if nargout>1
            u(i,:) = uk(:).';
        end
    end
end