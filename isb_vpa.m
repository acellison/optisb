function f = isb_vpa(x,Pfree,Pdep,Vfree,Vdep)
%function f = isb_vpa(pol)
%function f = isb_vpa(x,nfs,nds)
%function f = isb_vpa(x,Pfree,Pdep,Vfree,Vdep)
%   Compute the ISB.
    if nargin==1
        Pfree = 0;
        Pdep = fliplr(x);
        Vfree = 0;
        Vdep = 1;
        x = 0;
    elseif nargin==3
        nfs = Pfree;
        nds = Pdep;
        P = gbspolys(max(max(nfs),max(nds)),class(vpa(0)));
        Pfree = P(nfs,:);
        Pdep = P(nds,:);
        Vfree = rot90(vandern(1./vpa(nfs).^2,length(nds)));
        Vdep = rot90(vander(1./vpa(nds).^2));
    end

    if isempty(x); x=0; end
    if isempty(Pfree)||isempty(Vfree); Pfree=0; Vfree=0; end
    
    x = x(:);
    if length(x) ~= size(Pfree,1)
        error('Size of x must match the free thread polynomial table height');
    end
    if length(x) ~= size(Vfree,2)
        error('Size of x must match the width of its Vandermonde matrix');
    end
    if size(Vdep,1) ~= size(Vdep,2)
        error('Dependent Vandermonde matrix must be square');
    end
    if size(Vdep,1) ~= size(Pdep,1)
        error('Height of dependent thread polynomial table must match Vandermonde matrix size');
    end
    
    f = 0;

    % Algorithm tolerances
    itol = 1e-6;    % roots considered purely real if |imag(root)|<itol
    utol = 1e-12;   % roots a,b considered identical if max(|a-b|)<utol

    % Compute the coefficients of the threads
    b = [1;zeros(size(Vdep,1)-1,1,class(vpa(0)))];
    c = Vdep\(b-Vfree*x);

    % Compute our stability polynomial r(y) = |R(iy)|^2-1
    R = x.'*Pfree+c.'*Pdep;
    [r,~,~] = make_stability_poly(R);
    r = fliplr(r);
    
    % Evaluate the polynomial to ensure we are either negative or zero and
    % concave down at y=0
    r0 = vpapolyval(r,0);
    if r0<0 || (r0==0 && vpaconcud(r,0)<0)
        % Construct the companion matrix of to our stability polynomial
        rcon = condition_poly(r);
        M = compan(rcon);

        % Compute the eigenvalues
        zers = eig(M);
        
        % Throw away complex roots, as these aren't an imaginary axis
        % crossing
        inds = abs(imag(zers))<itol;
        zers = real(zers(inds));

        % Remove negative roots- the polynomial is symmetric about the
        % origin
        inds = zers>=0;
        zers = zers(inds);

        % Sort the roots in ascending order
        zers = sort(zers);

        % Find the first root where the function changes sign
        i = 1;
        while i<=length(zers)
            mult = 1;
            for j=i+1:length(zers)
                if abs(zers(j)-zers(i)) <= utol*max(abs(zers(:)))
                    mult = mult+1;
                else
                    break
                end
            end
            if mod(mult,2)==1
                indx = i;
                break
            end
            i = i+mult;
        end
        f = zers(indx);
    end
end

function [r,p,q] = make_stability_poly(R)
    % r = |R(iy)|^2 - 1 = p(y)^2 + q(y)^2 - 1
    % R input in reverse polyval order
    % r,p,q returned in reverse polyval order
    n = length(R)-1;
    p = zeros(size(R),class(R));
    q = zeros(size(R),class(R));
    if mod(n,2) == 0
        p(1:2:n+1) = (-1).^(0:n/2).*R(1:2:n+1);
        q(2:2:n)   = (-1).^(0:n/2-1).*R(2:2:n);
    else
        p(1:2:n)   = (-1).^(0:(n-1)/2).*R(1:2:n);
        q(2:2:n+1) = (-1).^(0:(n-1)/2).*R(2:2:n+1);
    end
    r = vpaconv(p,p)+vpaconv(q,q);
    r(1) = r(1)-1;
end

function r = condition_poly(c)
    % Taken from roots.m, precondition the polynomial for eigenvalue
    % computation
    c = c(:).';

    inz = find(c);
    if isempty(inz)
        % All elements are zero
        return
    end

    % Strip leading zeros and throw away.  
    % Strip trailing zeros, but remember them as roots at zero.
    nnz = length(inz);
    c = c(inz(1):inz(nnz));

    % Prevent relatively small leading coefficients from introducing Inf
    % by removing them.
    d = c(2:end)./c(1);
    while any(isinf(d))
        c = c(2:end);
        d = c(2:end)./c(1);
    end

    r = [1 d];
end

