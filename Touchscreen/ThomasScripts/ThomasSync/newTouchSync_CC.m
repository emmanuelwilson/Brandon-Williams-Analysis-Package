function [] = newTouchSync_CC(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
for i = 1 : length(folders)
    folders{i}
    if ~isempty(folders{i})        
        d = dir(folders{i});
        fnames = {d.name};
        if isempty(find(strcmp(fnames,'msTouchSync_new.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp.dat',1),1))
            cd(folders{i});                                  %Change current folder
            try
                single_folder_synchronization_V4()                                                                         %Run analysis
            catch
                fprintf([folders{i},' Failed to analize'])
            end            
        end
    end
end
end