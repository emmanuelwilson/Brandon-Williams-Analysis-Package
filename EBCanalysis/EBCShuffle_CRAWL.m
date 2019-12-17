%% Crawl script for EBC shuffling
%Must have all wanted folders under your "Current Folder" 
clear all
name = 'EBCresultsEvenOddSplitParallelNEW_deconvolvedFIX';
folder = dir(pwd);                                                         %List contents of your current folder
oldCD = pwd;                                                               %Save path to directory
for i = 3 : length(folder)                                                 %Look through folder items
    if folder(i).isdir == 1                                                %Will proceed only if folder item is another folder/directory
        subdir = [pwd,'\',folder(i).name];                                 %file path
        subfolder = dir(subdir);                                           %List subfolder contents
        fnames = {subfolder.name};                                         %List subfolder item names
        if ~isempty(find(strcmp(fnames,'msDeconvolved.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp',1),1))&& ~isempty(find(strncmp(fnames,'HeadTrackingData.mat',1),1))&& ~isempty(find(strncmp(fnames,'frameMap.mat',1),1)) && ~isempty(find(strncmp(fnames,name,1),1))
            cd([pwd,'/',folder(i).name]);                                  %Change current folder
            p = pwd;            
            f2 = dir(pwd);
            for ii = 3 : length(f2)                                                 %Look through folder items                
                subdir2 = [pwd,'\',f2(ii).name];                                 %file path
                subfolder2 = dir(subdir2);                                           %List subfolder contents
                fnames2 = {subfolder2.name};                                         %List subfolder item names                
                if f2(ii).isdir == 1 && strcmp(f2(ii).name, name) && ~isempty(find(strncmp(fnames2,'EBCstats.mat',1),1))
                    load('msDeconvolved.mat');
                    load('HeadTrackingData.mat');
                    load('frameMap.mat');
                    cd([pwd,'/',f2(ii).name]);                                  %Change current folder
                    load('EBCstats.mat')                    
                    msEgoCentricRateMapSplitEvenOddParallelSHUFFLE_V3(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimY, 1, out.QP)
                    break
                end
            end            
            cd(oldCD);                                                     %Return
        end
    end
end