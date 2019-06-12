function y = gbspolys(me,klass)
    if nargin==1
        klass = 'double';
    end
    me = ceil(double(me));
    n = me+1;
    y = zeros(n,n+1,klass);
    y(1,1:2) = 1;        % Build up the polynomials corresponding to equation
    y(2,1:3) = [1 2 2];  % (4.2) in [FZL].
    for k = 3:n
        y(k,1:k) = y(k-2,1:k);
        y(k,2:k+1) = y(k,2:k+1)+2*y(k-1,1:k);
    end
    for k = 2:2:n-1
        y(k,:) = 0.25*(y(k-1,:)+2*y(k,:)+y(k+1,:));
    end
    % Divide entries % Procedure here differs somewhat from [FZL].
    for k = 2:1:n
        y(k,:) = y(k,:)./(k.^(0:n));
    end
    % Remove entries on odd lines (for readability of y-array
    y(1:2:n,:) = 0;
end
