function []  = SIT_Start_Frame_Hab_Crawl(p)

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
        if ~isempty(find(strncmp(fnames,'behavCamCat.avi',15),1)) && ~isempty(find(strncmp(fnames,'frameMap.mat',10),1)) && isempty(find(strncmp(fnames,'ObjectStatsHab.mat',17),1))
            cd(folders{i});                                  %Change current folder
            load('frameMap.mat')
            load('badframes.mat')                      
            load('ObjectStats.mat')
            load('HeadTrackingData.mat')
            vid = VideoReader('behavCam1.avi');            
            frame = read(vid,t);
            h = figure;
            imshow(frame)
            hold on
            plot(SINKdata(frameMap(t:t+9000),1),SINKdata(frameMap(t:t+9000),2))
            [ObjectPos1, ObjectPos2] = getpts(h);
            Object1hab = cat(2,ObjectPos1,ObjectPos2);
            close(h);
            
            ObjectStatsHab.CentroidOb1 = Object1hab;
            save('ObjectStatsHab.mat','ObjectStatsHab')
            
            close all
            
        elseif ~isempty(find(strncmp(fnames,'behavCamCat.avi',15),1)) && ~isempty(find(strncmp(fnames,'frameMap.mat',10),1)) && ~isempty(find(strncmp(fnames,'ObjectStatsHab.mat',17),1))
%             cd(folders{i});                                  %Change current folder
%             load('ObjectStatsHab.mat')
%             load('SITstartFrameHab.mat')
%             load('HeadTrackingData.mat')
%             load('frameMap.mat')
%             h = figure;
%             plot(SINKdata(frameMap(startframe:end),1),SINKdata(frameMap(startframe:end),2))
%             [ObjectPos1, ObjectPos2] = getpts(h);
%             Object1hab = cat(2,ObjectPos1,ObjectPos2);
%             close(h);
%             
%             ObjectStatsHab.CentroidOb1 = Object1hab;
%             save('SITstartFrameHab.mat','startframe','ObjectStatsHab')
        end
    end
end