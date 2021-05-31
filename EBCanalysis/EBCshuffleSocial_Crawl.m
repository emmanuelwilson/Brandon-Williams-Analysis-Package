%% Sweetmilk/social target shuffling
function []  = EBCshuffleSocial_Crawl(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name};
        if ~isempty(find(strcmp(fnames,'ms.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp',1),1)) && ~isempty(find(strncmp(fnames,'EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3',1),1)) && isempty(find(strncmp(fnames,'ObjectShuffledStats.mat',15),1))
            cd(folders{i});                                  %Change current folder
            try
                load('EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3.mat/EBCstats.mat')
                load('ms.mat')
                load('HeadTrackingData.mat')
                load('frameMap.mat')
                msEgoCentricRateMapSplitEvenOddParallelSweetMilkSHUFFLE_V4_2(ms,HDdeg,SINKdata, frameMap, out.dimX, out.dimY, 1, out.QPO,out.QPW, 1,3);
            catch
                display([folders{i},' Failed to analize'])
            end
        end
    end
end
end