%% Fix indices from out of FOV cells and find/fill missing cell indices

function map = FindMissingCells2020(map, OOF,alignment)
ind = [];
count = 1;
if isempty(OOF)
    combs = find(~cellfun(@isempty,alignment.alignmentMap(1,:)));
    if ~isempty(combs)
        [~,~,NotmissingInd] = intersect(map(:,1),alignment.alignmentMap{1,combs(1)}(:,1));
        missingInd = find(alignment.alignmentMap{1,combs(1)}(:,1));
        missingInd = setdiff(missingInd,NotmissingInd);
        missingCellind = sort(alignment.alignmentMap{1,combs(1)}(missingInd,1),'descend');
        addon = zeros(length(missingCellind),length(map(1,:)));
        addon(:,1) = missingCellind;
        map = cat(1,map,addon);
        for i = 1 : length(combs)
            [~,~,NotmissingInd] = intersect(map(:,combs(i)),alignment.alignmentMap{1,combs(i)}(:,2));
            missingInd = find(alignment.alignmentMap{1,combs(i)}(:,2));
            missingInd = setdiff(missingInd,NotmissingInd);
            missingCellind = sort(alignment.alignmentMap{1,combs(i)}(missingInd,2),'descend');
            addon = zeros(length(missingCellind),length(map(1,:)));
            addon(:,combs(i)) = missingCellind;
            map = cat(1,map,addon);
        end
    end
else
    %update cell map for any out of FOV cells
    for i = 1 : length(OOF.exclude_all)
        valsind = find(map(:,i));
        realind = find(OOF.exclude_all{i});
        for j = 1 : length(valsind)
            map(valsind(j),i) = realind(map(valsind(j),i));
        end
    end
    
    %finds missing cells and adds index
    cellnum = 0;
    for i = 1 : length(map(1,:))
        missing = OOF.exclude_all{i};
        missing(map(find(map(:,i)),i)) = 0;
        missingCellind = find(missing);
        if ~isempty(missing)
            missingCellind = sort(missingCellind,'descend');
            addon = zeros(length(missingCellind),length(map(1,:)));
            addon(:,i) = missingCellind;
            map = cat(1,map,addon);
        end
    end
end
end