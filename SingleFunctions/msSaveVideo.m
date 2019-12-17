function  msSaveVideo(ms,frameLimit, downSample, columnCorrect, align, dFF, outFile)

b = load('frameMap_video.mat');
syncedFrames = b.frameMap;

%MSSAVEVIDEO Summary of this function goes here
%   Detailed explanation goes here
    mindFF = 0.1;
    maxdFF = 0.7;
    if isempty(frameLimit)
        frameLimit = [1 ms.numFrames];
    end
    
    writerObj = VideoWriter([outFile  '.avi'],'Grayscale AVI');
    open(writerObj);
    
    %for frameNum=frameLimit(1):downSample:frameLimit(2)
    for frameNum=1:length(syncedFrames)
        %frame = double(ms.mask).*msReadFrame(ms,frameNum,columnCorrect, align, dFF);
        frame = double(ms.mask).*msReadFrame(ms,syncedFrames(frameNum),columnCorrect, align, dFF);
        if dFF
            frame = 255*(frame-mindFF)/maxdFF;   
        else
            frame = 255*(frame-min(ms.minFluorescence))/max(ms.maxFluorescence);
        end
        
        writeVideo(writerObj,uint8(frame));
    
    end
    
    close(writerObj);
end

