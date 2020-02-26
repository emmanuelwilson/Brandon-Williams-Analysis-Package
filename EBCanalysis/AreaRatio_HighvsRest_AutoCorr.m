%% Takes ratio from Area under the curve from max peak vs everything else

function [ratios] = AreaRatio_HighvsRest_AutoCorr(corrangles,window)
maxval = [];
minval = [];
ratios = NaN(length(corrangles(:,1,1)),length(corrangles(1,:,1)));
y = [];
x = [];
% window = 8;
for i = 1 : length(corrangles(:,1,1))
    if any(any(~cellfun(@isempty,corrangles(i,:,:))))
        corrrots = permute(corrangles(i,:,:),[3 2 1]);        
        [y,x] = find(permute(~cellfun(@isempty,corrangles(i,:,:)),[3 2 1]));        
        for j = 1 : length(y)
            maxpeak = find(corrrots{y(j),x(j)} == max(corrrots{y(j),x(j)}));
            minamp = min((corrrots{y(j),x(j)}));
            if maxpeak < window+1
                topadd = window +1 - maxpeak;
                botadd = window - topadd;
                maxval = sum(corrrots{y(j),x(j)}([maxpeak-botadd:maxpeak+window, end - topadd : end]))-((window*2+1)*minamp);
                mincorr = corrrots{y(j),x(j)};
                mincorr([maxpeak-botadd:maxpeak+window, end - topadd : end]) = [];
                minval = sum(mincorr)-(length(mincorr)*minamp);
            elseif maxpeak > 60 - window
                botadd = maxpeak + window - 60;
                topadd = window-botadd;
                maxval = sum(corrrots{y(j),x(j)}([maxpeak-window:maxpeak+topadd, 1 : botadd]))-((window*2+1)*minamp);
                mincorr = corrrots{y(j),x(j)};
                mincorr([maxpeak-window:maxpeak+topadd, 1 : botadd]) = [];
                minval = sum(mincorr)-(length(mincorr)*minamp);
            else
                maxval = sum(corrrots{y(j),x(j)}(maxpeak-window:maxpeak+window))-((window*2+1)*minamp);
                mincorr = corrrots{y(j),x(j)};
                mincorr(maxpeak-window:maxpeak+window) = [];
                minval = sum(mincorr)-(length(mincorr)*minamp);
            end
            ratios(i,x(j)) = maxval/minval;
        end
    end
end

end