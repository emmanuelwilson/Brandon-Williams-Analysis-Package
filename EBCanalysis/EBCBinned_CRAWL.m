%% Crawl script for EBC shuffling
%Must have all wanted folders under your "Current Folder" 
clear all
name = 'EBCevenOddSplitParallelDeconvolved_CircSHUFFLEDonut3_100itt.mat';
namenogo = 'EBCresultsEvenOddSplitParallel_Deconvolved_Donut3_Bin_D1.2A6';
namenogo2 = 'EBCresultsEvenOddSplitParallel_Deconvolved_Donut3_Bin_D1A6';

folder = dir(pwd);                                                         %List contents of your current folder
oldCD = pwd;                                                               %Save path to directory
for i = 3 : length(folder)                                                 %Look through folder items
    if folder(i).isdir == 1                                                %Will proceed only if folder item is another folder/directory
        subdir = [pwd,'\',folder(i).name];                                 %file path
        subfolder = dir(subdir);                                           %List subfolder contents
        fnames = {subfolder.name};                                         %List subfolder item names
        if ~isempty(find(strcmp(fnames,'msDeconvolved.mat'),1)) && ~isempty(find(strncmp(fnames,'HeadTrackingData.mat',1),1)) && ~isempty(find(strncmp(fnames,'frameMap.mat',1),1)) && ~isempty(find(strncmp(fnames,name,60),1)) && isempty(find(strncmp(fnames,namenogo,59),1))&& isempty(find(strncmp(fnames,namenogo2,55),1))
            cd([pwd,'/',folder(i).name]);  
            load('msDeconvolved.mat');
            load('HeadTrackingData.mat');
            load('frameMap.mat');
            load(name)
            cd(subdir);
            if mod(out.dimX,1.2) == 0 && mod(out.dimY,1.2) == 0
                msEgoCentricRateMapSplitEvenOddParallel_V3_BINTEST(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimY, 1, out.QP,1:length(ms.FiltTraces(1,:)),1.2,6)
            else
                msEgoCentricRateMapSplitEvenOddParallel_V3_unevenBINs(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimY, 1, out.QP,1:length(ms.FiltTraces(1,:)),1.2,6)
            end
            cd(oldCD);                                                              %Return
        end
    end
end