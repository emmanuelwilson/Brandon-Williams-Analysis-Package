%% Crawl script for ms structure creation
%Must have all wanted folders and subfolders to run CNMFE on in the given path

function [] = msRun2020_oldSoft_CRAWL_V2(p)
paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
for i = 1 : length(folders)
%     folders{i}
    if ~isempty(folders{i})        
        d = dir(folders{i});
        fnames = {d.name};
        if isempty(find(strcmp(fnames,'ms.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp',1),1))
            cd(folders{i});                                  %Change current folder
            try
                msRun2018                                                      %Run analysis
            catch
                display([folders{i},' Failed to analize'])
            end            
        end
    end
end
end