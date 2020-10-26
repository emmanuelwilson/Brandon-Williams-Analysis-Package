%% Crawl script for ms structure creation
%Must have all wanted folders and subfolders to run CNMFE on in the given path

function [] = msRun2020_newSoft_CRAWL_V2_ZH(p)
paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
for i = 1 : length(folders)
    if ~isempty(folders{i})        
        d = dir(folders{i});
        fnames = {d.name};
        if isempty(find(strcmp(fnames,'ms.mat'),6)) && ~isempty(find(strncmp(fnames,'timeStamp',9),1)) && ~isempty(find(contains(folders{i},'Miniscope'),1))
            cd(folders{i});                                  %Change current folder
            msRun2020_newSoft_ZH(pwd) 
            %try
            %    msRun2020_newSoft(pwd)                                                      %Run analysis
            %catch
            %    display([folders{i},' Failed to analize'])
            %end            
        end
    end
end
end