%% Fix indices from out of FOV cells and find/fill missing cell indices

function map = FindMissingCells2020_Auto(map,sessions)
%finds missing cells and adds index
for i = 1 : length(map(1,:))
    missing = zeros(sessions{i}.numNeurons,1);
%     cellid = [1:sessions{i}.numNeurons];
    sind = find(map(:,i));
    if length(sind) < sessions{i}.numNeurons
        for c = 1 : sessions{i}.numNeurons
            miss = 0;
            for j = 1 : length(sind)
                if map(sind(j),i) == c
                    miss = 0;
                    break
                else 
                    miss = 1;
                end
            end
            missing(c) = miss;
        end        
    end  
    missingCellind = find(missing);
    if ~isempty(missing)
        missingCellind = sort(missingCellind,'descend');
        addon = zeros(length(missingCellind),length(map(1,:)));
        addon(:,i) = missingCellind;
        map = cat(1,map,addon);
    end
end
end