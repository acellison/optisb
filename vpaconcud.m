function ud = vpaconcud(p,x)
% function ud = vpaconcud(p,x)
% Check whether a polynomial is concave up or down at the point x.
% Returns +1 if concave up, -1 if concave down, and 0 otherwise.
    ud = 0;
    while length(p)>=3
        pp = vpapolyder(p);
        p  = vpapolyder(pp);
        if vpapolyval(pp,x)==0
            px = vpapolyval(p,x);
            if px<0
                ud = -1;
                break;
            elseif px>0
                ud = 1;
                break;
            end
        end
    end        
end
