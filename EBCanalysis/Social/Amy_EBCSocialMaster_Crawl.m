%% Sweetmilk/social target shuffling
function []  = Amy_EBCSocialMaster_Crawl(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
mouse = [];
frameMap1 = [];
frameMap2 = [];
HD1 = [];
HD2 = [];
tracking1 = [];
tracking2 = [];
ms = [];
out = [];
badframe1 = [];
badframe2 = [];
t1t2 = 0;

for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name};
        if (~isempty(find(strcmp(fnames,'ms.mat'),6)) && ~isempty(find(strncmp(fnames,'timestamp',9),1)) && ~isempty(find(strncmp(fnames,'EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3.mat',61),1)) && ~isempty(find(strncmp(fnames,'EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3_rightOb',63),1))) || (~isempty(find(strcmp(fnames,'ms.mat'),6)) && ~isempty(find(strcmp(fnames,'EBCstats.mat'),8)) && ~isempty(find(strcmp(fnames,'PassedVals.mat'),10)))
            cd(folders{i});                                  %Change current folder
            subfolders = strsplit(folders{i},'\');
            for j = length(subfolders) : -1 : 1
                try
                    if ~isempty(find(contains(subfolders{j}, 'camkii-'),1)) && isempty(find(contains(subfolders{j}, 'trial'),1)) && isempty(find(contains(subfolders{j}, 'habituation'),1)) || ~isempty(find(contains(subfolders{j}, 'Trial1Trial2ObjectData_D1A3_bothOb'),1))
                        newName = subfolders{j};
                        if isempty(find(contains(newName, 'camkii-'),1)) && ~isempty(find(contains(subfolders{j}, 'Trial1Trial2ObjectData_D1A3_bothOb'),1))
                            t1t2 = 1;
                            for s = length(subfolders) : -1 : 1
                                if ~isempty(find(contains(subfolders{s}, 'camkii-'),1))
                                    newName = subfolders{s};
                                end
                            end
                        end
                        if isempty(find(strcmp(newName, mouse),1))
                            mouse = newName;
                            day = subfolders{j-1};
                            frameMap1 = [];
                            frameMap2 = [];
                            HD1 = [];
                            HD2 = [];
                            tracking1 = [];
                            tracking2 = [];
                            ms = [];
                            out = [];
                            badframe1 = [];
                            badframe2 = [];                            
                        end
                        if t1t2
                            load('ms.mat')
                            load('EBCstats.mat')
                            t1t2 = 0;
                        end
                        break
                    elseif ~isempty(find(contains(subfolders{j}, 'trial'),1))
                        if ~isempty(find(contains(subfolders{j}, 'trial1'),1))
                            load('frameMap.mat')
                            load('HeadTrackingData.mat')
                            load('badframes.mat')
                            
                            frameMap1 = frameMap;
                            tracking1 = SINKdata;
                            HD1 = HDdeg;
                            badframe1 = t;
                        elseif ~isempty(find(contains(subfolders{j}, 'trial2'),1))
                            load('frameMap.mat')
                            load('HeadTrackingData.mat')
                            load('badframes.mat')
                            
                            frameMap2 = frameMap;
                            tracking2 = SINKdata;
                            HD2 = HDdeg;
                            badframe2 = t;
                        end
%                     elseif ~isempty(find(contains(subfolders{j}, 'Trial1Trial2ObjectData_D1A3_bothOb'),1))                        
                    end
                end
            end
        end
        if ~isempty(frameMap1) && ~isempty(frameMap2) && ~isempty(HD1) && ~isempty(HD2) && ~isempty(tracking1) && ~isempty(tracking2) && ~isempty(ms) && ~isempty(out)
            try
                cd('../../')
                msEgoCentricRateMapSplitEvenOddParallelSweetMilkParallel_V2(ms,HD1,HD2,tracking1,tracking2, frameMap1,frameMap2, out.dimX, out.dimY, 1, out.QPOL,out.QPOR, out.QPW,1,3,badframe1,badframe2)
                close all
                frameMap1 = [];
                frameMap2 = [];
                HD1 = [];
                HD2 = [];
                tracking1 = [];
                tracking2 = [];
                ms = [];
                out = [];
                badframe1 = [];
                badframe2 = [];
            catch
                display([folders{i},' Failed to analize'])
            end
        end
    end
end