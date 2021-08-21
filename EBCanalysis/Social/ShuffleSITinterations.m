
function out = ShuffleSITinterations(traces,enterzones,exitzones,shuffleNum)

out = zeros(length(traces(:,1)),shuffleNum,length(enterzones));
interactionframes = exitzones - enterzones; 
indshortint = [];
for i = 1 : shuffleNum
    shuffledtraces = CShuffle(traces');
    for j = 1: length(traces(:,1))
        for s = 1 : length(enterzones)
            if interactionframes(s) >= 60
                out(j,i,s) = mean(shuffledtraces(enterzones(s):exitzones(s),j))/interactionframes(s);
            else
                indshortint = cat(1,indshortint,s);
            end
        end
    end
end
indshortint = sort(unique(indshortint),'descend');
out(:,:,indshortint) = [];
end