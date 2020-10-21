%% Concactenate videos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT: 
%       -folderpath: Input the folder location with the videos you wish to
%       combine
%       -filePrefix: Will combine videos with said prefix
%   OUTPUT:
%       Will save videos in the folder path under the names of
%       [filePrefix]Cat.avi.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = catVids(folderpath, filePrefix)

folder = dir(folderpath);

%Count # of frames
aviFiles = dir([folderpath,'\', '*.avi']);
vidcount = 1;
numFrames = 0;
bigvid = [];
cmap = gray(256);
%Take all frames and make into one
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
%Write video
v = VideoWriter([filePrefix,'Cat.avi'],'Uncompressed AVI');
% v.Colormap = bigvid(1).colormap;
open(v)
for j = 1 : length(bigvid(:,1))
    writeVideo(v,bigvid(j).cdata)
end
close(v)
end