function result = ratapprox(x,delta,maxdigits)
    maxden = 10^maxdigits;
    s = sign(x);
    x = abs(x);
    scale = 2^floor(log2(x));
    x = x/scale;

    result = sym([]);
    for b=1:maxden
        alo =  ceil(b*(x*(1-delta)));
        ahi = floor(b*(x*(1+delta)));
        if alo<=ahi
            r = s*(sym(alo:ahi)/sym(b));
            result = [result r];
        end
    end
    result = unique(sort(result));
    result = scale*result;
end
