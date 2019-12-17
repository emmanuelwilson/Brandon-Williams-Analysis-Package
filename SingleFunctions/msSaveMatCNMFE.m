function  msSaveMatCNMFE(name,frameLimit, scale)
%MSSAVEVIDEO Summary of this function goes here
%   Detailed explanation goes here

vidObj = dir(name).name;
v = VideoReader(vidObj);
if isempty(frameLimit)
    frameLimit = [1 v.numFrames];
end

Y=[];
for frameNum=frameLimit(1):1:frameLimit(2)
    frame = msReadFrame(v,frameNum,0, 0, 0);
    if scale ~= 1
        frame = imresize(frame,scale);
    end
    
    Y(:,:,frameNum)=frame;    
    
    if (mod(frameNum,1000)==0)
        display(['Creating video, on frame: ' num2str(frameNum)])
    end
       
end

Ysiz = size(Y);
savefast('msvideo.mat','Y', 'Ysiz');

end

