%% Fix indices from out of FOV cells and find/fill missing cell indices

function map = FindMissingCells2020_Auto(map,OOF,alignment,combs)
ind = [];
count = 1;
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