%% Will create multiple videos with different
function [] = VidSplitBenchmark(folderpath)
format long
folder = dir(folderpath);

%User promtp for # of subvideos to make
% prompt = 'How many splits do you want?';
% splitnum = input(prompt,'s');
% [splitnum,isnumber] = str2num(splitnum);
% while ~isnumber
%     if isinteger(splitnum)
%         isnumber = true;
%     else
%         splitnum = input(prompt,'s');
%         [splitnum,isnumber] = str2num(splitnum);
%     end
% end

%Count # of frames
aviFiles = dir([folderpath,'\', '*.avi']);
filePrefix = 'msCam';
vidcount = 1;
numFrames = 0;
bigvid = [];
cmap = gray(256);
for i=1:length(aviFiles)
    endIndex = strfind(aviFiles(i).name,'.avi');
    if (~isempty(strfind(aviFiles(i).name,filePrefix)))
        vidObj{vidcount} = VideoReader([folderpath,'\', filePrefix, num2str(vidcount), '.avi']);
        numFrames = numFrames + vidObj{vidcount}.NumberOfFrames;
        b = VideoReader([folderpath,'\', filePrefix, num2str(vidcount), '.avi']);
        if isempty(bigvid)
            bigvid = im2frame(readFrame(b),cmap);
            for j = 2 : vidObj{vidcount}.NumberOfFrames
                bigvid(j,:) = im2frame(readFrame(b),cmap);
            end
        else
            for j = length(bigvid(:,1))+1 : (length(bigvid(:,1))) + vidObj{vidcount}.NumberOfFrames
                bigvid(j,:) = im2frame(readFrame(b),cmap);
            end
        end
        vidcount = vidcount + 1;
    end
end

%Look at main FOV and create subFOVs
v = VideoReader([folderpath,'\msCam1.avi']);
refFrame = readFrame(v(1));
if mean(mean(refFrame))<3
    refFrame = readFrame(v);
end
figure
imshow(refFrame)

dims = [round(length(refFrame(1,:))*0.6),round(length(refFrame(:,1))*0.6)];

[x,y] = getpts;
xmin = zeros(length(x),1);
ymin = zeros(length(y),1);
imsec = zeros(length(x),4);
for i = 1 : length(x)
    xmin(i) = x(i) - round(dims(1)/2);
    ymin(i) = y(i) - round(dims(2)/2);
    [~,imsec(i,:)] = imcrop(refFrame, [xmin(i), ymin(i), dims(1), dims(2)]);
end

sL = round(numFrames/length(x));
if sL * length(x) > numFrames
    spill = true;
else
    spill = false;
end
for i = 1 : length(x)
    v = VideoWriter([filePrefix,'Split',num2str(i),'.avi'],'grayscale AVI');   
    open(v)
    if  spill && i == length(x)
        over = sL * length(x) - numFrames;
        for j = 1 : sL-over
            writeVideo(v,bigvid(int32(j+(sL*(i-1)))).cdata(int16(imsec(i,2)):int16(imsec(i,2))+int16(round(imsec(i,4))),int16(imsec(i,1)):int16(imsec(i,1))+int16(round(imsec(i,3)))))
        end
    else
        for j = 1 : sL
            writeVideo(v,bigvid(int32(j+(sL*(i-1)))).cdata(int16(imsec(i,2)):int16(imsec(i,2))+int16(round(imsec(i,4))),int16(imsec(i,1)):int16(imsec(i,1))+int16(round(imsec(i,3)))))
        end
    end
    close(v)
end
end