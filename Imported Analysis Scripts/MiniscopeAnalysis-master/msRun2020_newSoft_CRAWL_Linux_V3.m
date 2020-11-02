%% Crawl script for ms structure creation
%Must have all wanted folders and subfolders to run CNMFE on in the given path
%This version is made for Compute Canada analysis on a UNIX system and will
%convert all analysis calcium videos to GREY before running analysis

function [] = msRun2020_newSoft_CRAWL_Linux_V3(p)
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
            cmdStr = sprintf('%s' , '/lustre03/project/6049321/m3group/Wilson/Brandon-Williams-Analysis-Package/SLURM/convert_grey_V2_updated.sh');            
            system(cmdStr);
            try            
                msRun2020_newSoft(pwd)                                                      %Run analysis
            catch
                display([folders{i},' Failed to analize'])
            end            
        end
    end
end
end