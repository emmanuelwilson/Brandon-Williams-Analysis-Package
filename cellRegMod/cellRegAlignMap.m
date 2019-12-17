%% Recovers cell registration and sorts matrices in cell array
function alignmentmap = cellRegAlignMap(folderpath)

folder = dir(folderpath);
n = length(folder(:,1))-2;
tn = roots([1 -1 -n*2]);
trialNum = tn(1);
alignmentmap = cell(int32(trialNum));

for i = 1 : trialNum-1
    for j = i+1 : trialNum
        ftemp = dir([folderpath,'\',num2str(i),'_',num2str(j)]);
        for s = 1 : length(ftemp)
            if contains(ftemp(s).name,'cellRegistered')
                reg = load([folderpath,'\',num2str(i),'_',num2str(j),'\',ftemp(s).name]);
                alignmentmap{i,j} = reg.cell_registered_struct.cell_to_index_map;
            end
        end
    end
end