function [RandPeaks] = ShufflePeaks(firing)%
%%Takes firing trace matrix and shuffles peak locations through time
%
%INPUT: -n by m matrix "firing", where n is the number of frames and m is 
%       the cell number.
%
%OUTPUT:-RandPeaks, n by m matrix (as defined with firing) containing all
%       peaks with amplitudes above 0.0001 in randomly selected points in time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Author: Emmanuel Wilson

firing(firing<0.0001) = 0;                                                  %Sets all values below 0.0001 to zero
Shuff = zeros(length(firing(:,1)),length(firing(1,:)));                     %Shuffled firing matrix
Peaks = zeros(500,length(firing(1,:)));                                     %Peak info matrix,values represents the number of frames the peak contains
indPeaks = Peaks;                                                           %Peak Index matrix,values are the index of the first frame in the peak

for cellNum = 1: length(firing(1,:))                                        %Going through each cell
    count = 0;                                                              %Counter for number of frames in each peak
    for i = 1 : length(firing(:,1))                                         %Going through each frame
        if firing(i,cellNum)>0 && count == 0                                %if first frame in peak
            count = count + 1;                                              %update counter
            indPeaks(length(find(indPeaks(:,cellNum)))+1,cellNum) = i;      %Assign frame# to peak index
        elseif firing(i,cellNum)>0                                          %if within the peak
            count = count + 1;                                              %update counter
        elseif count>0                                                      %if first frame after peak
            Peaks(length(find(Peaks(:,cellNum)))+1,cellNum)=count;          %record total counts for that peak
            count = 0;                                                      %update counter
        end
    end
    for j=1: length(find(Peaks(:,cellNum)))                                 %Going through each peak
        r = randi(length(firing(:,1)));                                     %random frame value between 1 and input
        while (r+Peaks(j,cellNum))> length(firing(:,1)) || any(Shuff(r:(r+Peaks(j,cellNum)-1),cellNum)) %will keep looping unless the peak can fit in its new location
            r = randi(length(firing(:,1)));                                 %new random frame value
        end
        Shuff(r:(r+Peaks(j,cellNum)-1),cellNum) = firing(indPeaks(j,cellNum):(indPeaks(j,cellNum)+Peaks(j,cellNum)-1),cellNum); %Places peak at random frame r.
    end
end
RandPeaks = Shuff;                                                          %Output
end