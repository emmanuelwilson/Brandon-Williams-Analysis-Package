%% Crawl script for deconvolved ms structure creation

clear all

folder = dir(pwd);
oldCD = pwd;
for i = 3 : length(folder)
    if folder(i).isdir == 1
        subdir = [pwd,'\',folder(i).name];
        subfolder = dir(subdir);
        fnames = {subfolder.name};
        if ~isempty(find(strcmp(fnames,'msTouchSync.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp.dat',1),1)) && ~isempty(find(strcmp(fnames,'msDeconvolved.mat'),1))
            cd([pwd,'/',folder(i).name]);
            load('msDeconvolved.mat')
            load('msTouchSync.mat')
            
            TUNLtaskFigureGenerator_Optomized_V2_Figcentric
            
            cd(oldCD);
        end
    end
end