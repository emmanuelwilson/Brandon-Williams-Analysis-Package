function [] = newTouchSync_CRAWL()

folder = dir(pwd);
oldCD = pwd;
for i = 3 : length(folder)
    if folder(i).isdir == 1
        subdir = [pwd,'\',folder(i).name];
        subfolder = dir(subdir);
        fnames = {subfolder.name};
        if isempty(find(strcmp(fnames,'msTouchSync_new.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp.dat',1),1))
            cd([pwd,'/',folder(i).name]);
            csvFiles = dir('*.csv');
            for j = 1 : length(csvFiles)
                if ~isempty(find(strncmp(fnames,csvFiles(j).name,1),1))                    
                    try
                        single_folder_synchronization_V4()
                    end
                end
            end
            cd(oldCD)
        end
    end
end