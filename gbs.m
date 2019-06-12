function uk = gbs(f,h,t0,u0,n)
% function yk = gbs(f,h,t0,u0,n)
% Run N steps of the GBS algorithm on function f with step size h, starting
% at initial condition (t0,u0).  n must be even so that even-order accuracy
% requirements are met.
    if mod(n,2)~=0
        error(['Number of steps is not a multiple of two.  '...
               'Even-order accuracy is not guaranteed']);
    end
    u0 = u0(:).';
    u  = zeros(n+2,length(u0),class(u0));
    t  = t0;
    u(1,:) = u0;
    u(2,:) = u(1,:)+h*f(t,u(1,:));                % Forward Euler Step
    for i=3:n+2
        t = t+h;
        u(i,:) = u(i-2,:)+2*h*f(t,u(i-1,:));      % Leap Frog Iteration
    end
    uk = .25*(u(end-2,:)+2*u(end-1,:)+u(end,:));  % Smoothing Step
end
