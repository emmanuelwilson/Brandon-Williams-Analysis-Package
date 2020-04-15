function [pos uninterp] = trackPosition_1LED(doPath)
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

    frame  = 1;                 %Frame # of current video
    CurrentFrame = frame;
    pos = nan(2, totalFrames);
    maxVal = nan(1, totalFrames);
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
                thisFrame = read(videoObject, frame);
                thisFrame = thisFrame(:,:,1);
                maxVal(CurrentFrame) = nanmax(thisFrame(:));
                if nanmax(thisFrame(:)) < 150
                    frame = frame + 1;
                    CurrentFrame = CurrentFrame + 1;
                    continue
                end
                [x y] = find(thisFrame == nanmax(thisFrame(:)));
                pos(:,CurrentFrame) = [nanmedian(x) nanmedian(y)]';
                frame = frame + 1;
                CurrentFrame = CurrentFrame + 1;
            end
            toc
            fprintf('\b')
        end
        frame = 1;
    end
    uninterp = pos;
    pos = interpNaNs(pos')';
    
    figure(1)
    set(gcf,'position',[50 50 1.*nanmax(pos')-nanmin(pos')])
    plot(pos(1,:),pos(2,:),'color',[0.2 0.2 0.2])
    axis equal
    axis off
    qi = find(ismember(doPath,'/'));
    doPath(qi) = '_';
    outP = ['Plots/Paths/' doPath(qi(2)+1:qi(3)-1) '/' doPath(qi(2)+1:end)];
    saveFig(gcf,outP,'tiff')
    close all
end