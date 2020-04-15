function [pos uninterp] = trackPosition_2LED(doPath)
folder = dir(doPath);     %Looks at current file location

    totalFrames = 0;        %Total number of frames observed (including previous video itterations)
    TimestampMin = 10000000;%Arbitrary VERY large number for min value storage, must be larger than datenum stamp in folder
    Timestamp = 0;          %TimeStamp record for indexing
    o = NaN(1,length(folder)); %index movie order

    % find avi and dat files
    aviFiles = dir([doPath '\behavCam*.avi']);
    filePrefix = 'behavCam';
    
    clipNum = nan(1,length(aviFiles));
    for i = 1:length(aviFiles)
        clipNum(i) = str2num(aviFiles(i).name(9:end-4));
    end
    [a b] = sort(clipNum);
    aviFiles = aviFiles(b);
    clipNum = clipNum(b);

    %Find total number of frames and organize vidoes by their timestamps 
    for ind_vid = 1: length(folder)                         %Going through every object in the folder
        if (~isempty(strfind(folder(ind_vid).name,filePrefix)))
            movieFullFileName =  folder(ind_vid).name;      %Movie file name
            video = VideoReader([doPath '/' movieFullFileName]);   %Read video file              
            totalFrames = totalFrames + video.NumberOfFrames;    %Determine how many frames in file
            for t = 1 : length(folder)
                if (~isempty(strfind(folder(t).name,filePrefix)) && folder(t).datenum < TimestampMin && folder(t).datenum > Timestamp)
                    TimestampMin = folder(t).datenum;       %Save as minimum value
                    o(1,ind_vid) = t;                       %Save Index #
                end
            end
            Timestamp = TimestampMin;                       %Set minimum time boundary
            TimestampMin = 10*TimestampMin;                 %Reset minimum value
        end               
    end
    [mx my] = meshgrid(1:video.Width,1:video.Height);
    
    frame  = 1;                 %Frame # of current video
    CurrentFrame = frame;
    pos = nan(2, totalFrames, 2);
    maxVal = nan(1, totalFrames);
    sensitivityThresh = [50 50];
    for it = 1 : length(o)                 %Going through every object in the folder
        ind_vid = o(1,it);
        if ~isnan(ind_vid) && (folder(ind_vid).isdir == 0)                     %ignore if file is not a video
            movieFullFileName =  folder(ind_vid).name;      %Movie file name
            videoObject = VideoReader([doPath '/' movieFullFileName]);   %Read video file      
            numberOfFrames = videoObject.NumberOfFrames;    %Determine how many frames in current video file   
            
            
            v = VideoReader([doPath '/' movieFullFileName]);             %Read video file a gain ( hasFrame() function won't work otherwise)

            fprintf(['\n\t\tTracking Frames (Clip ' num2str(it-2) '):  '])
            
            tic
            while hasFrame(v) && frame<= numberOfFrames
                ogFrame = read(videoObject, frame);
                for i = 1:2
                    if i == 1
                        thisFrame = ogFrame(:,:,i) - nanmax(ogFrame(:,:,[1:i-1 i+1:end]),[],3);
                        rFrame = thisFrame;
                    else
                        thisFrame = ogFrame(:,:,i) - rFrame;
                    end
                    maxVal(CurrentFrame) = nanmax(thisFrame(:));
                    if nanmax(thisFrame(:)) < sensitivityThresh(i)
                        continue
                    end
                    df = double(thisFrame);
                    df(df<(nanmax(df(:))./2)) = 0;
                    df = df./nansum(df(:));
                    x = nansum(df(:).*mx(:));
                    y = nansum(df(:).*my(:));
%                     [x y] = find(thisFrame == nanmax(thisFrame(:)));
                    pos(:,CurrentFrame,i) = [nanmedian(x) nanmedian(y)]';
                end
                frame = frame + 1;
                CurrentFrame = CurrentFrame + 1;
            end
            toc
            fprintf('\b')
        end
        frame = 1;
    end
    pos(:,pos(1,:,1)==0 & pos(2,:,1)==0,1) = nan;
    pos(:,pos(1,:,2)==0 & pos(2,:,2)==0,2) = nan;
    uninterp = pos;
    pos(:,:,1) = interpNaNs(pos(:,:,1)')';
    pos(:,:,2) = interpNaNs(pos(:,:,2)')';
    
    figure(1)
%     set(gcf,'position',[50 50 1.*nanmax(pos(:,:,1)')-nanmin(pos(:,:,1)')])
    plot(pos(1,:,1),pos(2,:,1),'color',[0.8 0.5 0.5])
    hold on
    plot(pos(1,:,2),pos(2,:,2),'color',[0.5 0.8 0.5])
    axis equal
    axis off
    qi = find(ismember(doPath,'/'));
    doPath(qi) = '_';
    outP = ['Plots/Paths/' doPath(qi(2)+1:qi(3)-1) '/' doPath(qi(2)+1:end)];
    saveFig(gcf,outP,'tiff')
    close all
end