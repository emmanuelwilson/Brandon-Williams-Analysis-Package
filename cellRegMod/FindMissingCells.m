%% Fix indices from out of FOV cells and find/fill missing cell indices

function map = FindMissingCells(map, folder)
resultsfolder = dir([folder,'\Results']);
ind = [];
count = 1;
%update cell map for any out of FOV cells
for i = 1 : length(resultsfolder)
    if resultsfolder(i).isdir == 0 && contains(resultsfolder(i).name,'cellmap')
        session = str2num(resultsfolder(i).name(1+7:end-4));
        load([resultsfolder(i).folder,'\',resultsfolder(i).name])
        vals = find(map(:,session));
        for j = 1 : length(vals)
            map(vals(j),session) = cellmap(map(vals(j),session));
        end
    end
end

%finds missing cells and adds index
folderdir = dir(folder);
cellnum = 0;
for i = 1 : length(map(1,:))
    for j = 1 : length(folderdir)
        if (folderdir(j).isdir == 0) && (contains(folderdir(j).name,['ms',num2str(i),'.mat']))
            load([folder,'/ms',num2str(i),'.mat'])
            cellnum = 1:ms.numNeurons;
            if length(find(map(:,i))) < ms.numNeurons
                vals = map(find(map(:,i)),i);
                vals = sort(vals,'descend');
                cellnum(vals) = [];
                addon = zeros(length(cellnum),length(map(1,:)));
                addon(:,i) = cellnum;
                map = cat(1,map,addon);
                continue
            end
        end
    end
end

% %Fix indices for out of FOV cells
% if ~isempty(ind)
%     for i = 1 : length(ind)
%         session = str2num(resultsfolder(ind(i)).name(1+7:end-4));
%         load([resultsfolder(ind(i)).folder,'\',resultsfolder(ind(i)).name])
%         vals = find(map(:,session));
%         for j = 1 : length(vals)
%             map(vals(j),session) = cellmap(map(vals(j),session));
%         end
%     end
% end



%Add missing cell indices
end