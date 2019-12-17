function output = MovieSplitScreen(dirName)

filePrefix = 'behavCam';
% load('msTouchSync.mat')

% for i = 1 : length(events(:,1))
%     frameMap(i,1) = str2num(char(events(i,1)));
% end

MAXFRAMESPERFILE = 1000; %This is set in the miniscope control software

% find avi and dat files
aviFiles = dir([dirName '\*.avi']);
datFiles = dir([dirName '\*.dat']);
brain = dir('msvideo.avi');
folder = dir(dirName);

ms.numFiles = 0;        %Number of relevant .avi files in the folder
ms.numFrames = 0;       %Number of frames within said videos
ms.vidNum = [];         %Video index
ms.frameNum = [];       %Frame number index
ms.maxFramesPerFile = MAXFRAMESPERFILE; %finds the maximum number of frames contained in a single throughout all videos
ms.dffframe = [];

%find the total number of relevant video files
for i=1:length(aviFiles)
    endIndex = strfind(aviFiles(i).name,'.avi');        %find the name of current .avi file
    if (~isempty(strfind(aviFiles(i).name,filePrefix)))
        ms.numFiles = max([ms.numFiles str2double(aviFiles(i).name((length(filePrefix)+1):endIndex))]);     % +1 count for relevant .avi files
    end
end

o = NaN(1,length(folder));
for i = 1:length(folder)
    if ~folder(i).isdir && ~isempty(strfind(folder(i).name,filePrefix))
        o(i) = str2num(folder(i).name(9:end-4));
    end
end
[a b] = sort(o);
folder = folder(b);
o = o(b);

BrainRead = VideoReader([brain.folder '\' brain.name]);
i = 1;
while hasFrame(BrainRead)
    video(:,:,i) = readFrame(BrainRead);    
    i = i+1;
end
resizeFlag = 0;
for vidCount = 1 : ms.numFiles
    j = o(1,vidCount);
    BehavRead = VideoReader(folder(j).name);     
    i = 1+1000*(vidCount-1);
    while hasFrame(BehavRead)        
        Beh(:,:,i) = rgb2gray(readFrame(BehavRead));
        if ~(size(Beh(:,:,i)) == size(video(:,:,1)))
            BehVid(:,:,i) = imresize(Beh(:,:,i),[length(video(:,1,1)) length(video(1,:,1))]);
            resizeFlag = 1;
        end
        i = i+1;
    end
end
if resizeFlag == 0
    BehVid = Beh;
end


SplitFrame = cat(2,video,BehVid(:,:,frameMap));

v = VideoWriter('SplitVidSpeedUp.avi');
open(v)
for i = 1 : length(frameMap(:,1))
    writeVideo(v,SplitFrame(:,:,i))
    i = i+9;
end
close(v)

% 
% open(v)
% for i = 1 : totalFrames
%     BrainRead = VideoReader([brain.folder '\' brain.name]);
%     Imbrain = read(BrainRead,i);
%     if frameMap(i)>1000*vidCount
%         vidCount = vidCount+1;
%         j = o(vidCount,1);        
%     end
%     BehavRead = videoreader(folder(j).name);
%     ImBeh = readFrame(frameMap(i)-(1000*(vidCount-1)));
%     writeVideo(v,cat(2,Imbrain,ImBeh));
% end
% close(v)

end