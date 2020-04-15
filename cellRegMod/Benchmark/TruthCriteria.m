%Compare truth with results
totbadreg = 0;
totgoodreg = 0;
totmissedreg = 0;
perstats = zeros(length(tmap(:,1)), 3);
missedind = [];
count = 0;
ind = [];
for i = 1 : length(tmap(:,1))
    rows = zeros(length(tmap(1,:)),1);
    
    start = length(find(tmap(i,:)));
    if start < length(tmap(1,:)) && isempty(find(tmap(i,:),1))
        c = find(tmap(i,:),1);
    else
        c = 1;
    end
    for j = c: length(tmap(1,:))
        if j > 1 && tmap(i,j) == 0 %rows(j-1) == find(Singlemap(:,j) == tmap(i,j))
            rows(j) = rows(j-1);
            m = 1;
            while rows(j-m) < 0
                m = m +1;
                rows(j) = rows(j-m);
            end
            if length(Singlemap(rows(j),j) == tmap(i,j)) < 2 && Singlemap(rows(j),j) == tmap(i,j)
                
            else
                rows(j) = -rows(j-1);
            end
        else
            rows(j) = find(Singlemap(:,j) == tmap(i,j));
        end
    end
    
    numgvals = 0;
    numbvals = 0;
    missedvals = 0;
    sess2 = false;
    if i == 14
    end
    rowsame = mode(rows);
    if length(find(rows == rowsame)) == length(tmap(1,:))
        perstats(i,1) = 1;
        perstats(i,2) = 0;
        perstats(i,3) = 0;
        numgvals = length(tmap(1,:));
    else
        urow = unique(rows);
        missedvals = length(urow) - 1;
        for j = 1 :length(urow)
            if find(urow(j)< 0)
                numbvals = numbvals + length(find(rows< 0));
                missedvals = missedvals -1;
            elseif find(urow(j) == 0)
                if rowsame > 0
                    rcomp = rowsame;
                else
                    rcomp = rows(find(rows> 0,1));
                end
                if rcomp > 0 && isempty(find(Singlemap(rcomp,1:length(find(rows == 0))),1))
                    numgvals = numgvals + length(find(rows == 0));
                    missedvals = missedvals - 1;
                else
                    numbvals = numbvals + length(find(Singlemap(rcomp,1:length(find(rows == 0)))));
                    numgvals = length(find(rows == 0)) - length(find(Singlemap(rcomp,1:length(find(rows == 0)))))+1;
                    missedvals = missedvals - 1;
                end
            elseif length(find(rows == urow(j))) == 1 && length(find(abs(rows) == urow(j))) > 1
%                 missedvals = missedvals + 1;
            elseif length(find(rows == urow(j))) == 1
                if length(find(Singlemap(urow(j),:))) == 1
%                     missedvals = missedvals + 1;
                else
%                     numbvals = numbvals + 1;
                end
            else
                sess = length(find(rows == urow(j)));
                numgvals = numgvals + sess;
                if sess2
                    numgvals = numgvals -1;
                end
                sess2 = true;
                if length(find(Singlemap(urow(j),:))) > sess
                    CorrectwIncorrect(i,j) = i;
                end
            end
        end
        perstats(i,1) = numgvals/length(tmap(1,:));
        perstats(i,2) = numbvals/length(tmap(1,:));
        perstats(i,3) = missedvals/length(tmap(1,:));
    end
    if (numbvals + numgvals + missedvals) > 12 || (numbvals + numgvals + missedvals) < 12
        count = count + 1;
        ind(count) = i;
    end
    totbadreg = totbadreg + numbvals;
    totgoodreg = totgoodreg + numgvals;
    totmissedreg = totmissedreg + missedvals;
end