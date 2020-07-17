%% Crawl script for ms structure creation
%Must have all wanted folders and subfolders to run CNMFE on in the given path

function [] = msRun2020_newSoft_CRAWL(p)

paths = genpath(startPath);
folders = strsplit(paths,';')';

for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name};
        if isempty(find(strcmp(fnames,'ms.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp',1),1))
            cd([p,'/',folder(i).name]);                                  %Change current folder
            try
                msRun2020_newSoft(pwd)                                                      %Run analysis
            catch
                display([p,'/',folder(i).name,' Failed to analize'])
            end
            cd(oldCD);
        end
    end
end
end