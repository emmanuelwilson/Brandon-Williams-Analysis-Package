%% Crawling script for EBC HeadTracking 

clear all
folder = dir(pwd);
oldCD = pwd;
for i = 3 : length(folder)
    if folder(i).isdir == 1
        subdir = [pwd,'\', folder(i).name];
        subfolder = dir(subdir);
        fnames = {subfolder.name};
        if isempty(find(strcmp(fnames,'HeadTrackingData.mat'),1))
            cd(subdir);
            HDTrackingE
            cd(oldCD);
        end
    end
end