function z = vpaconv(x,y)
    if ~isa(x,'sym'), x = cast(x,'sym'); end
    if ~isa(y,'sym'), y = cast(y,'sym'); end
    N = length(x);
    M = length(y);
    z = zeros(1,N+M-1,'sym');

    f = [x zeros(1,M,class(x))];
    g = [y zeros(1,N,class(y))];

    for i=1:length(z)
        for j=1:length(x)
            if i-j+1>0
                z(i) = z(i) + f(j)*g(i-j+1);
            end
        end
    end
end
