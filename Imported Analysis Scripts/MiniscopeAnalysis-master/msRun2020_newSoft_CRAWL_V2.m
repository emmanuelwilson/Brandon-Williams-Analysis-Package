%% Crawl script for ms structure creation
%Must have all wanted folders and subfolders to run CNMFE on in the given path

function [] = msRun2020_newSoft_CRAWL_V2(p)
paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
folders
for i = 1 : length(folders)    
    if ~isempty(folders{i})        
        folders{i}
        d = dir(folders{i});
        fnames = {d.name};
        if isempty(find(strcmp(fnames,'ms.mat'),1)) && ~isempty(find(contains(folders{i},'Miniscope'),1)) %&& ~isempty(find(strncmp(fnames,'timeStamps.csv',9),1)) && ~isempty(find(contains(folders{i},'Miniscope'),1))
            cd(folders{i});                                  %Change current folder
            try
                msRun2020_newSoft(pwd)                                                      %Run analysis
            catch
                display([folders{i},' Failed to analize'])
            end            
        end
    end
end
end