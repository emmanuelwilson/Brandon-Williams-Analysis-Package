function [nmap,avg_psame] = ElimConflict2020(map,alignmentMap,probmap)

count = 0;
ecount = 0;
elim = [];
avg_psame = zeros(length(map(:,1)),1);

%average probability of same registration across all sessions
for i = 1 : length(map(:,1))
    comb = map(i,:);
    loc = find(map(i,:));
    if length(loc) > 1
        for k = 1 : length(loc)-1
            fsess = true;
            q = k+1;
            while fsess
                tempRegMap = alignmentMap{loc(k), loc(q)};
                temp1 = find(tempRegMap(:,1) == comb(loc(k)));
                if tempRegMap(temp1,2) == comb(loc(q))
                    tempprob = probmap{loc(k), loc(q)};
                    val = (tempprob(temp1,1));
                    count = count+ 1;
                    avg_psame(i) = avg_psame(i) + val;
                    fsess = false;
                end
                q = q +1;
                if q > length(loc)
                    fsess = false;
                end
            end
            
        end
        avg_psame(i) = avg_psame(i)/count;
        count = 0;
    else
        avg_psame(i) = 0;
    end
end

%find conflicting cell registration from last to first session
for i = length(map(1,:)) : -1 : 1
    totcell = length(unique(find(map(:,i))));
    %For each cell
    for j = 1: totcell
        pairs = find(map(:,i) == j);
        if length(pairs)>1
            pSame = avg_psame(pairs);
            [~, I] = max(pSame);            
            pairs(I) = [];
            if isempty(elim)
                elim = pairs;
            else
                elim(end+1:end+length(pairs)) = pairs;
            end
            pairs = [];
        end        
    end
end

elim = unique(elim);
nmap = map;
avg_psame(elim,:) = [];
nmap(elim,:) = [];
end