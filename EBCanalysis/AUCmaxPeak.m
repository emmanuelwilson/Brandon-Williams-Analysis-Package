%% Takes Area under the curve from a user defined distance from the max peak

function [AUC] = AUCmaxPeak(corrangles,window)
maxval = [];
minval = [];
AUC = NaN(size(corrangles));
y = [];
x = [];

for i = 1 : length(corrangles(:,1,1))
    if any(any(~cellfun(@isempty,corrangles(i,:,:))))
        corrrots = permute(corrangles(i,:,:),[3 2 1]);
        [y,x] = find(permute(~cellfun(@isempty,corrangles(i,:,:)),[3 2 1]));
        for j = 1 : length(y)
            maxpeak = find(corrrots{y(j),x(j)} == max(corrrots{y(j),x(j)}));
            if maxpeak < window+1
                topadd = window +1 - maxpeak;
                botadd = window - topadd;
                maxval = sum(corrrots{y(j),x(j)}([maxpeak-botadd:maxpeak+window, end - topadd : end]));                
            elseif maxpeak > 52
                botadd = maxpeak + window - 60;
                topadd = window-botadd;
                maxval = sum(corrrots{y(j),x(j)}([maxpeak-window:maxpeak+topadd, 1 : botadd]));               
            else
                maxval = sum(corrrots{y(j),x(j)}(maxpeak-window:maxpeak+window));                
            end
            AUC(i,y(j),x(j)) = maxval;
        end
    end
end

end