function [RandPeaks] = ShuffleSpikeTrain(firing)%, frameMap, SINKdata, HDdeg)
%%Takes firing trace matrix and shuffles peak locations through time
%
%INPUT: -n by m matrix spike train of neural activity, where n is the number
%       of frames and m is the cell number.
%
%OUTPUT:-RandPeaks, n by m matrix (as defined with firing) containing all
%       spikes in randomly selected points in time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Author: Emmanuel Wilson

Shuff = zeros(length(firing(:,1)),length(firing(1,:)));                     %Shuffled firing matrix

for cellNum = 1: length(firing(1,:))                                        %Going through each cell
    indPeaks = find(firing(:,1));                                           %Counter for number of frames in each peak
    for i = 1 : length(indPeaks(:,1))                                       %Going through each frame
        r = randi(length(firing(:,1)));                                     %random frame value between 1 and input
        while (r> length(firing(:,1)) || any(Shuff(r,cellNum)))             %will keep looping unless the peak can fit in its new location
            r = randi(length(firing(:,1)));                                 %new random frame value
        end
        Shuff(r,cellNum) = firing(indPeaks(i),cellNum);                     %Places peak at random frame r.
    end
end
RandPeaks = Shuff;                                                          %Output
end