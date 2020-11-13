%% Reads videos in folder and outputs average frame intensity

function [meanintensity] = MeanFrameIntensity(path)

oldcd = pwd;
cd(path)

avidir = dir([pwd '/*.avi']);

meanintensity = [];

for i = 1 : length(avidir)
    for j = 1 : length(avidir)
    if strcmp([num2str(i-1) '.avi'], avidir(j).name)
        vid = VideoReader(avidir(j).name);
        count = length(meanintensity) + 1;
        meanintensity = cat(1,meanintensity, zeros(vid.NumFrames,1));
        while hasFrame(vid)
            frame = readFrame(vid);
            frame = rgb2gray(frame);
            meanintensity(count) = (mean(mean(frame)));
            count = count +1;
        end
    end
    end
end
end