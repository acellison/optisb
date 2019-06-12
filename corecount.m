function argout = corecount(ns,kk)
% function argout = corecount(ns,kk)
%
% Compute the minimum number of cores required to run the given step
% counts, provided no core takes more than max(ns) function calls.
% If kk is provided, checks we can fit on at least kk cores.
    ns = ns(:);
    ns = sort(ns,'descend');
    nmax = ns(1);
    m = length(ns);
    argout = m;
    for k=1:m-1
        fits = trykcores(ns(2:end),k,nmax);
        if fits
            argout = k+1;
            break;
        end
    end
    if nargin>1
        if argout<=kk
            argout = 1;
        else
            argout = 0;
        end
    end
end

function fits = trykcores(ns,k,nmax)
    sums = zeros(k,1);
    used = zeros(length(ns),1);
    for i=1:k
        for j=1:length(ns)
            if ~used(j)
                tmp = sums(i)+ns(j);
                if tmp<=nmax
                    sums(i) = tmp;
                    used(j) = 1;
                end
            end
        end
    end
    if all(used==1)
        fits = 1;
    else
        fits = 0;
    end
end