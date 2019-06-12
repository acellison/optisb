function [uk,us] = rexstep(scheme,f,h,t0,u0,ns,cs,n)
% function uk = rexstep(scheme,f,h,t0,u0,ns,cs,n)
% function [uk,us] = rexstep(scheme,f,h,t0,u0,ns,cs,n)
% Run N richardson extrapolations, performing individual time steps with
% scheme.
% param scheme: callable as scheme(f,h/ns(i),t0,u0,ns(i))
% param f:  callable as f(xi,yi)
% param t0: initial time point
% param u0: initial function value
% param ns: vector of time steps performed in the extrapolated threads
% param cs: weights of each extrapolated thread
% param n:  number of extrapolated time steps to perform.  Defaults to 1
% return uk: array of extrapolated outputs at times t0+h*i, i=1:n
    if nargin==7
        n = 1; 
    end

    u0 = u0(:).';
    ns = ns(:);
    cs = cs(:);
    if length(ns) ~= length(cs)
        error('Step size counts must match number of extrapolation coefficients');
    end

    klass = class(u0);
    if nargout>1
        us = zeros(n,length(u0),klass);
    end
    
    t    = t0;
    uk   = u0;
    uext = zeros(length(ns),length(u0),klass);

    % Sort the weights in increasing order to reduce roundoff error in
    % summing threads with drastically disparate weights.
    schemes = cell(length(cs),1);
    for i=1:length(cs)
        schemes{i} = scheme;
    end
    hs = h./ns;
    for i=1:n
        % Run the time stepper routines.  This can be performed in parallel
        for j=1:length(ns)
            uext(j,:) = schemes{j}(f,hs(j),t,uk,double(ns(j)));
        end
            
        % Combine each of the threads to produce the output sample at time t+h
        uk = cast(cs.'*uext,klass);

        % Store the intermediate result if desired
        if nargout>1
            us(i,:) = uk;
        end
        
        % Step the time index
        t = t+h;
    end
end

