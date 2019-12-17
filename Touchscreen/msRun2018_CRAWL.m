%% Crawl script for ms structure creation
%Must have all wanted folders under your "Current Folder" 
clear all

folder = dir(pwd);                                                         %List contents of your current folder
oldCD = pwd;                                                               %Save path to directory
for i = 3 : length(folder)                                                 %Look through folder items
    if folder(i).isdir == 1                                                %Will proceed only if folder item is another folder/directory
        subdir = [pwd,'\',folder(i).name];                                 %file path
        subfolder = dir(subdir);                                           %List subfolder contents
        fnames = {subfolder.name};                                         %List subfolder item names
        if isempty(find(strcmp(fnames,'ms.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp',1),1)) %Will proceed if no structure is present and timestamp file being present
            cd([pwd,'/',folder(i).name]);                                  %Change current folder
            msRun2018                                                      %Run analysis 
            cd(oldCD);                                                     %Return
        end
    end
end