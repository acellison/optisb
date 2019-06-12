function cells = sortbyfield(cells,field)
    result = zeros(size(cells),class(cells{1}.(field)));
    for i=1:length(cells)
        result(i) = cells{i}.(field);
    end
    [~,idx] = sort(result);
    cells = {cells{idx}};
end