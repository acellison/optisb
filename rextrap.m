function [pol, weights] = rextrap(polys,ns,rhs)
% function [pol, weights] = rextrap(polys,ns)
% function [pol, weights] = rextrap(polys,ns,rhs)
% Richardson Extrapolate the polynomials with corresponding step count
% sizes.  Set up the Vandermonde system and solve to eliminate successive
% terms in the error expansion, so that we eliminate the lowest
% (length(ns))-order terms in the error expansion.  All rows in polys
% must be generated using the same time-stepping method to guarantee each
% coefficient ak for error term ak*h^k in the asymptotic expansions match
% so that these can be eliminated.  The returned polynomial has the same
% order as the largest order of the input polynomials.
    if size(polys,1)~=length(ns)
        error('Number of polynomials must match the step count length');
    end
    if nargin>2
        b = rhs;
    else
        b = [1;zeros(size(polys,1)-1,1)];
    end
    V = rot90(vander(1./ns));
    c = V\b;
    pol = (polys'*c)';
    if nargout>1
        weights = c;
    end
end