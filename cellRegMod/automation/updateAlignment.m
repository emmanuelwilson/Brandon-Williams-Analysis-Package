%% Will reasign the cell indecies to take into account the cells which were removed

function alignmentnew = updateAlignment(alignment,OOF,combs)
for i = 1 : length(combs)
    if ~isempty(alignment{i})
        inds = find(alignment{i}(:,2));
        vals = alignment{i}(inds,2);
        rvals = OOF.cellmap{i}(vals);
        for j = 1 : length(find(alignment{i}(:,2)))
            alignmentnew{combs(i,1),combs(i,2)} = alignment{i};
            alignmentnew{combs(i,1),combs(i,2)}(inds(j)) = rvals(j);
        end
    else
        alignmentnew{combs(i,1),combs(i,2)} = [];
    end
end


end