function p = vpapolyval(pol,z)
    p = pol(1)*ones(size(z),class(pol));
    for i=2:length(pol)
        p = p.*z+pol(i);
    end
end
