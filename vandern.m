function V = vandern(v,n)
    v = v(:);
    if isempty(v)
        V = reshape(v, 0, n);
    else
        V = repmat(v, 1, n);
        V(:, n) = 1;
        V = cumprod(V, 2, 'reverse');
    end
end