%% Reverse Crawl script for EBC shuffling
%Must have all wanted folders under your "Current Folder" 
clear all

folder = dir(pwd);                                                         %List contents of your current folder
oldCD = pwd;                                                               %Save path to directory
for i = length(folder) : -1 : 3                                                 %Look through folder items
    if folder(i).isdir == 1                                                %Will proceed only if folder item is another folder/directory
        subdir = [pwd,'\',folder(i).name];                                 %file path
        subfolder = dir(subdir);                                           %List subfolder contents
        fnames = {subfolder.name};                                         %List subfolder item names
        if ~isempty(find(strcmp(fnames,'msDeconvolved.mat'),1)) && ~isempty(find(strncmp(fnames,'frameMap',1),1)) && ~isempty(find(strncmp(fnames,'HeadTrackingData',1),1)) && isempty(find(strncmp(fnames,'EBCevenOddSplitParallelDeconvolved_CircSHUFFLEDonut3_100itt',59),1))
            if ~isempty(find(strncmp(fnames,'EBCresultsEvenOddSplitParallelNEW_Deconvolved_Donut3',1),1))                
                cd([pwd,'/',folder(i).name]);                                  %Change current folder
                subdir2 = [pwd,'\EBCresultsEvenOddSplitParallel_Deconvolved_Donut3'];                                 %file path
                subfolder2 = dir(subdir2);                                           %List subfolder contents
                fnames2 = {subfolder2.name};                                         %List subfolder item names
                if ~isempty(find(strcmp(fnames2,'EBCstats.mat'),1)) %&& isempty(find(strcmp(fnames2,'EBCevenOddSplitParallelDeconvolved_CircSHUFFLE.mat'),1))
                    load('EBCresultsEvenOddSplitParallel_Deconvolved_Donut3/EBCstats.mat')
                    load('msDeconvolved.mat')
                    load('frameMap.mat')
                    load('HeadTrackingData.mat')
                    msEgoCentricRateMapSplitEvenOddParallelSHUFFLE_V3(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimY, 1, out.QP)                                                      %Run analysis
                end
                cd(oldCD);            
            elseif ~isempty(find(strncmp(fnames,'EBCresultsEvenOddSplitParallelNEW_deconvolvedFIX',1),1))
                cd([pwd,'/',folder(i).name]);                                  %Change current folder
                subdir2 = [pwd,'\EBCresultsEvenOddSplitParallelNEW_deconvolvedFIX'];                                 %file path
                subfolder2 = dir(subdir2);                                           %List subfolder contents
                fnames2 = {subfolder2.name};                                         %List subfolder item names
                if ~isempty(find(strcmp(fnames2,'EBCstats.mat'),1)) %&& isempty(find(strcmp(fnames2,'EBCevenOddSplitParallelDeconvolved_CircSHUFFLE.mat'),1))
                    load('EBCresultsEvenOddSplitParallelNEW_deconvolvedFIX/EBCstats.mat')
                    if length(fieldnames(out)) == 16
                        load('msDeconvolved.mat')
                        load('frameMap.mat')
                        load('HeadTrackingData.mat')
                        msEgoCentricRateMapSplitEvenOddParallelSHUFFLE_V3(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimY, 1, out.QP)                                                      %Run analysis
                    else
                        load('msDeconvolved.mat')
                        load('frameMap.mat')
                        load('HeadTrackingData.mat')
                        msEgoCentricRateMapSplitEvenOddParallelSHUFFLE_V3(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimX, 1, out.QP)                                                      %Run analysis
                    end
                end
                cd(oldCD);
            elseif ~isempty(find(strncmp(fnames,'EBCresultsEvenOddSplitParallelNEW_deconvolved',1),1))
                cd([pwd,'/',folder(i).name]);                                  %Change current folder
                subdir2 = [pwd,'\EBCresultsEvenOddSplitParallelNEW_deconvolved'];                                 %file path
                subfolder2 = dir(subdir2);                                           %List subfolder contents
                fnames2 = {subfolder2.name};                                         %List subfolder item names
                if ~isempty(find(strcmp(fnames2,'EBCstats.mat'),1)) %&& isempty(find(strcmp(fnames2,'EBCevenOddSplitParallelDeconvolved_CircSHUFFLE.mat'),1))
                    load('EBCresultsEvenOddSplitParallelNEW_deconvolved/EBCstats.mat')
                    if length(fieldnames(out)) == 16
                        load('msDeconvolved.mat')
                        load('frameMap.mat')
                        load('HeadTrackingData.mat')
                        msEgoCentricRateMapSplitEvenOddParallelSHUFFLE_V3(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimY, 1, out.QP)                                                      %Run analysis
                    else
                        load('msDeconvolved.mat')
                        load('frameMap.mat')
                        load('HeadTrackingData.mat')
                        msEgoCentricRateMapSplitEvenOddParallelSHUFFLE_V3(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimX, 1, out.QP)                                                      %Run analysis
                    end
                end
                cd(oldCD);
            elseif ~isempty(find(strncmp(fnames,'EBCresultsEvenOddSplitParallelNEWFix',1),1))
                cd([pwd,'/',folder(i).name]);                                  %Change current folder
                subdir2 = [pwd,'\EBCresultsEvenOddSplitParallelNEWFix'];                                 %file path
                subfolder2 = dir(subdir2);                                           %List subfolder contents
                fnames2 = {subfolder2.name};                                         %List subfolder item names
                if ~isempty(find(strcmp(fnames2,'EBCstats.mat'),1)) %&& isempty(find(strcmp(fnames2,'EBCevenOddSplitParallelDeconvolved_CircSHUFFLE.mat'),1))
                    load('EBCresultsEvenOddSplitParallelNEWFix/EBCstats.mat')
                    if length(fieldnames(out)) == 16
                        load('msDeconvolved.mat')
                        load('frameMap.mat')
                        load('HeadTrackingData.mat')
                        msEgoCentricRateMapSplitEvenOddParallelSHUFFLE_V3(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimY, 1, out.QP)                                                      %Run analysis
                    else
                        load('msDeconvolved.mat')
                        load('frameMap.mat')
                        load('HeadTrackingData.mat')
                        msEgoCentricRateMapSplitEvenOddParallelSHUFFLE_V3(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimX, 1, out.QP)                                                      %Run analysis
                    end
                end
                cd(oldCD);
            end
                                                   %Return
        end
    end
end