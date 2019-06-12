function isbn = isbn_est(pol)    
    fnz = find(pol~=0,1);
    order = length(pol)-fnz;

    fn = @(y) abs(polyval(fliplr(pol),order*complex(0,y)))-1;
    isbn = 0;
    if fn(0.1)*fn(1.0)<0; isbn = fzero(fn,[0.1,1.0]); end
end
