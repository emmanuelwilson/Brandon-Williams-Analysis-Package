%% Creates correlation map at the cell level with

function cellcorrmat = SingleCellCorrMap(map,alignmentMap, probmap)
cellcorrmat = zeros(length(map(1,:)),length(map(1,:)),length(map(:,1)));
bmap = map;
bmap(bmap>0) = 1;
bmap = sum(bmap,2);

for i = 1 : length(bmap)
    if bmap(i) > 1
        regloc = find(map(i,:));
        reg = map(i,regloc);
        d = zeros(length(map(1,:)),1);
        d(regloc) = 1;
        cellcorrmat(:,:,i) = diag(d);
        for j = 1 : length(regloc)-1
            for k = j + 1 : length(regloc)
                tempRegMap = alignmentMap{regloc(j), regloc(k)};
                temp1 = find(tempRegMap(:,1) == reg(j));
                if tempRegMap(temp1,2) == reg(k)
                    tempprob = probmap{regloc(j), regloc(k)};
                    val = (tempprob{temp1,1}(2));
                    cellcorrmat(regloc(j),regloc(k),i) =val;
                    cellcorrmat(regloc(k),regloc(j),i) =val;
                end
            end
            
        end
        
    end
    
end
end