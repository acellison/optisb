function maketight(ax,alpha)
    if nargin<2
        alpha=1;
    end
    inset = get(ax, 'TightInset');
    set(ax, 'Position', [inset(1), inset(2), alpha*(1-(inset(1)+inset(3))), alpha*(1-inset(2)-inset(4))])
end
