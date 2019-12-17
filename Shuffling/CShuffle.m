function [RandPeaks] = CShuffle(firing)

Shuff = zeros(length(firing(:,1)),length(firing(1,:)));                     %Shuffled firing matrix
minShift = 900;                                                             %minimum shift 900 frames: 30seconds

for cellNum = 1: length(firing(1,:))                                        %Going through each cell   
    
    r = randi([minShift length(firing(:,1))-minShift]);                                     %random frame value between 1 and input
    
    Shuff(:,cellNum) = circshift(firing(:,cellNum),r);                     %Shifts trace at random frame r.
end
RandPeaks = Shuff;                                                          %Output
end