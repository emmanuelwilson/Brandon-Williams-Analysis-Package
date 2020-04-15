%% Will compare the truth registration map with the test registration map and give stats on how good it did

totbadreg = 0;
totgoodreg = 0;
totmissedreg = 0;

for i = 1 : length(tmap(:,1))
    if ~isempty(find(tmap(i,:)))
        if ~isempty(find(ismember(Singlemap, tmap(i,:),'rows')))
            totgoodreg = totgoodreg + length(tmap(1,:));
        else
            rows = zeros(length(tmap(1,:)),1);
            for j = 1 : length(tmap(1,:))
                if tmap(i,j) > 0
                    rows(j) = find(Singlemap(:,j) == tmap(i,j));
                end
            end
            numgvals = 0;
            numbvals = 0;
            missedvals = 0;
            rowsame = mode(rows(rows>0));
            urow = unique(rows);
            
            if length(find(rows)) < length(rows)
                espace = find(rows==0);                
                tespace = find(tmap(i,:) == 0);
                ecommon = intersect(tespace, espace);
                
                numgvals = numgvals + length(ecommon);
                
                if length(tespace) > length(ecommon)
                    numbvals = numbvals + (length(tespace) - length(ecommon));
                end
                
            end
            urow = urow(urow>0);
            urow(find(urow == rowsame)) = [];
            if ~isempty(urow)
                for irow = 1 : length(urow)
                    if length(find(Singlemap(urow(irow),:))) == 1
                        missedvals = missedvals +1;
                    else
                        numbvals = numbvals + 1;
                    end
                end
            end
        end
        totbadreg = totbadreg + numbvals;
        totgoodreg = totgoodreg + numgvals;
        totmissedreg = totmissedreg + missedvals;
        numgvals = 0;
        numbvals = 0;
        missedvals = 0;
    end   
end