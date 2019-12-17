%% Moving average of values that are only backwards in time, ie not influenced by future events.

function [out] = BackMovingAvg(calcium)
windowSize = 30;                                                           %Window size in frames
out = zeros(size(calcium));
for i = 1 : length(calcium(:,1))
    for j = 1 : length(windowSize)
        if (i - windowSize) > 1
            out(i-j+1,:) = out(i-j+1,:) + mean(calcium(i-(windowSize-j+1):i-j+1,:),1);
        else
            sWindow = i - windowSize;
            if j >= sWindow
                break
            else
                out(i-j+1,:) = out(i-j+1,:) + mean(calcium(i-(sWindow-j+1):i-sWindow,:),1);
            end
        end
    end    
end
end