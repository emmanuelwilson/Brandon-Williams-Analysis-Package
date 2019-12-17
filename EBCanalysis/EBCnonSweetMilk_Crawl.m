%% Crawl script for EBC shuffling
%Must have all wanted folders under your "Current Folder" 
clear all
name = 'EBCevenOddSplitParallelDeconvolved_CircSHUFFLEDonut3SWEETMILK_D1A6.mat';
% namenogo = 'EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1d2A6.mat';
folder = dir(pwd);                                                         %List contents of your current folder
oldCD = pwd;                                                               %Save path to directory
for i = 3 : length(folder)                                                 %Look through folder items
    if folder(i).isdir == 1                                                %Will proceed only if folder item is another folder/directory
        subdir = [pwd,'\',folder(i).name];                                 %file path
        subfolder = dir(subdir);                                           %List subfolder contents
        fnames = {subfolder.name};                                         %List subfolder item names
        if ~isempty(find(strcmp(fnames,'msDeconvolved.mat'),1)) && ~isempty(find(strncmp(fnames,'HeadTrackingData.mat',1),1)) && ~isempty(find(strncmp(fnames,'frameMap.mat',1),1)) && ~isempty(find(strncmp(fnames,name,65),1)) %&& isempty(find(strncmp(fnames,namenogo,61),1))
            cd(subdir);
            load('msDeconvolved.mat');
            load('HeadTrackingData.mat');
            load('frameMap.mat');
            load(name)            
            msEgoCentricRateMapSplitEvenOddParallel_V3_BINTEST(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimY, 1, out.QPW,1:ms.numNeurons,1.2,6)
            cd(oldCD);                                                              %Return
        end
    end
end