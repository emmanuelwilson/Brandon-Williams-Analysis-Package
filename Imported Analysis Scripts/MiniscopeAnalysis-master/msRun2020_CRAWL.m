%% Crawl script for ms structure creation
%Must have all wanted folders and subfolders to run CNMFE on in the given path

function [] = msRun2020_CRAWL(p)

folder = dir(p);                                                         %List contents of your current folder
oldCD = pwd;                                                               %Save path to directory
for i = 3 : length(folder)                                                 %Look through folder items
    if folder(i).isdir == 1                                                %Will proceed only if folder item is another folder/directory
        subdir = [p,'/',folder(i).name];                                 %file path
        subfolder = dir(subdir);                                           %List subfolder contents
        fnames = {subfolder.name};                                         %List subfolder item names
        if isempty(find(strcmp(fnames,'ms.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp',1),1)) %Will proceed if no structure is present and timestamp file being present
            cd([p,'/',folder(i).name]);                                  %Change current folder
            try
                msRun2020(pwd)                                                      %Run analysis
            catch
                display([p,'/',folder(i).name,' Failed to analize'])
            end
            cd(oldCD);                                                     %Return
            
        else
            for j = 3: length(subfolder)
                if subfolder(j).isdir ==1
                    subdir2 = [p,'/',folder(i).name,'/',subfolder(j).name];                                 %file path
                    subfolder2 = dir(subdir2);                                           %List subfolder contents
                    fnames2 = {subfolder2.name};                                         %List subfolder item names
                    if isempty(find(strcmp(fnames2,'ms.mat'),1)) && ~isempty(find(strncmp(fnames2,'timestamp',1),1)) %Will proceed if no structure is present and timestamp file being present
                        cd([p,'/',folder(i).name,'/',subfolder(j).name]);                                  %Change current folder
                        try
                            msRun2020(pwd)                                                      %Run analysis
                        catch
                            display([p,'/',folder(i).name,'/',subfolder(j).name,' Failed to analize'])
                        end
                        msRun2020(pwd)                                                  %Run analysis
                        cd(oldCD);
                    end
                end
            end
        end    
    end
end
end