function vidObj = msDFFvidSave(vidObj,stepSize, frameLimits)
%   This function works to save a DF/F video from physiology recordings of
% calcium transience. The inputs are: (the ms file taken from msRun, the 
% step size (input 1 for a frame to frame video, 2 for double time, etc.),
% and the final input is a frame limit vector in case you don't want to
% film the entire physiology. If you do want to film the entire physiology
% then set frameLimits to [].
% Name your video file on line 33 and 34 or rename it after processing
if isempty(frameLimits) 
    frameLimits = [1 vidObj.numFrames];     %sets frame limits to the last frame
end

%smoothing kernal
hSmall = fspecial('average', 3);
hLarge = fspecial('average', 60);  % was 60

red = cat(3, ones(vidObj.alignedHeight(vidObj.selectedAlignment),vidObj.alignedWidth(vidObj.selectedAlignment)), ...
    zeros(vidObj.alignedHeight(vidObj.selectedAlignment),vidObj.alignedWidth(vidObj.selectedAlignment)), ...
    zeros(vidObj.alignedHeight(vidObj.selectedAlignment),vidObj.alignedWidth(vidObj.selectedAlignment)));

%used for overlay in display
green = cat(3, zeros(vidObj.alignedHeight(vidObj.selectedAlignment),vidObj.alignedWidth(vidObj.selectedAlignment)), ...
    ones(vidObj.alignedHeight(vidObj.selectedAlignment),vidObj.alignedWidth(vidObj.selectedAlignment)), ...
    zeros(vidObj.alignedHeight(vidObj.selectedAlignment),vidObj.alignedWidth(vidObj.selectedAlignment)));

% kernal used for eroding and dilating
se = strel('diamond',10);
se2 = strel('diamond',2);

frame = nan(vidObj.alignedHeight(vidObj.selectedAlignment),vidObj.alignedWidth(vidObj.selectedAlignment),stepSize);
vidObj.brightSpots = zeros(vidObj.alignedHeight(vidObj.selectedAlignment),vidObj.alignedWidth(vidObj.selectedAlignment)); %each cell (or pixel) contains the count of local maxima detected
vidObj.brightSpotTiming = sparse(vidObj.alignedHeight(vidObj.selectedAlignment)*vidObj.alignedWidth(vidObj.selectedAlignment),0); %holds the spatiotemporal information of detected bright spots
myVideo = VideoWriter('testy2_0.8.avi');
%uncompressedVideo = VideoWriter('testing4', 'Uncompressed AVI');
open(myVideo);

for startFrameNum=frameLimits(1):stepSize:min([frameLimits(2) vidObj.numFrames-stepSize]) %steps through data video
    count = 0;
    % Generates a max projection of the frames contained in the current
    % step. Applies a small spatial filter to the max projection to remove noise.
    % Applies a large spatial filter to get a measure of the background
    % activity.
    for frameNum = startFrameNum:(startFrameNum+stepSize) 
        count = count+1;
        frame(:,:,count) =filter2(hSmall,(double(vidObj.mask).*msReadFrame(vidObj,frameNum,true,true,true)));                        
    end
    frameMax = max(frame,[],3);
    frameBase = filter2(hLarge,frameMax);
    frameBase(frameBase<0) = 0;
    film = frameMax-frameBase;
    film(film<0) = 0;
    film(film>1) = 1;
    film = imadjust(film, [0 1], [0 1],0.8);                %attempt to brighten the video, comment out if necessary
    writeVideo(myVideo, film)
end
close(myVideo)
%save('testing')