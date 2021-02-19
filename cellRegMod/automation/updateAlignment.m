%% Will reasign the cell indecies to take into account the cells which were removed

function alignmentnew = updateAlignment(alignment,OOFm,OOFr,combs)
for i = 1 : length(combs)
    if ~isempty(alignment{i})
        %reference indexing
        indsr = find(alignment{i}(:,1));
        valsr = alignment{i}(indsr,1);
        rvalsr = OOFr.cellmap{i}(valsr);
        %moved indexing 
        indsm = find(alignment{i}(:,2));
        valsm = alignment{i}(indsm,2);
        rvalsm = OOFm.cellmap{i}(valsm);           
        
        alignmentnew{combs(i,1),combs(i,2)} = alignment{i};
        
        for j = 1 : length(find(alignment{i}(:,1)))            
            alignmentnew{combs(i,1),combs(i,2)}(indsr(j),1) = rvalsr(j);
        end
        for j = 1 : length(find(alignment{i}(:,2)))
            alignmentnew{combs(i,1),combs(i,2)}(indsm(j),2) = rvalsm(j);
        end
    else
        alignmentnew{combs(i,1),combs(i,2)} = [];
    end
end


end