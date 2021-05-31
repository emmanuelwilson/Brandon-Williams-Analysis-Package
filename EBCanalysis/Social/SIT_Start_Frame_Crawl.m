function []  = SIT_Start_Frame_Crawl(p)

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
        if ~isempty(find(strncmp(fnames,'behavCamCat.avi',15),1)) && ~isempty(find(strncmp(fnames,'frameMap.mat',10),1)) && isempty(find(strncmp(fnames,'SITstartFrame.mat',10),1))
            cd(folders{i});                                  %Change current folder
            load('frameMap.mat')
            lmov = length(frameMap);
            implay('behavCamCat.avi',100)
            startframe1 = input('Which frame does SIT start?');
            startframe = find(frameMap == startframe1);
            if isempty(startframe)
                startframe = find(frameMap == startframe1+1);
            end
            load('ObjectStats.mat')
            load('HeadTrackingData.mat')
            h = figure;
            plot(SINKdata(frameMap(startframe:end),1),SINKdata(frameMap(startframe:end),2))
            [ObjectPos1, ObjectPos2] = getpts(h);
            Object1 = cat(2,ObjectPos1,ObjectPos2);
            close(h);
            
            ObjectStats.CentroidOb1 = Object1;
            save('SITstartFrame.mat','startframe','ObjectStats')
            
            close all
            
        elseif ~isempty(find(strncmp(fnames,'behavCamCat.avi',15),1)) && ~isempty(find(strncmp(fnames,'frameMap.mat',10),1)) && ~isempty(find(strncmp(fnames,'SITstartFrame.mat',10),1))
            cd(folders{i});                                  %Change current folder
            load('ObjectStats.mat')
            load('SITstartFrame.mat')
            load('HeadTrackingData.mat')
            load('frameMap.mat')
            h = figure;
            plot(SINKdata(frameMap(startframe:end),1),SINKdata(frameMap(startframe:end),2))
            [ObjectPos1, ObjectPos2] = getpts(h);
            Object1 = cat(2,ObjectPos1,ObjectPos2);
            close(h);
            
            ObjectStats.CentroidOb1 = Object1;
            save('SITstartFrame.mat','startframe','ObjectStats')
        end
    end
end