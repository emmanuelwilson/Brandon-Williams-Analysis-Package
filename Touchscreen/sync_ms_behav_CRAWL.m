%% Crawl script for frameMap creation
clear all
folder = dir(pwd);
oldCD = pwd;
for i = 3 : length(folder)
    if folder(i).isdir == 1
        subfolder = dir([pwd,'/',folder(i).name]);
        fnames = {subfolder.name};
        if isempty(find(strncmp(fnames,'frameMap',1),1)) && ~isempty(find(strncmp(fnames,'timestamp',1),1))            
            cd([pwd,'/',folder(i).name]);
            sync_ms_behav
            cd(oldCD);
        end
    end
end