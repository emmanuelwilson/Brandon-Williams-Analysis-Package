%% will eliminate cell according to index at every 100 itterations.
function [newShuffle,Shuffle] = CellElimShuffle(Shuffle,elimInd)

for i = 1 : 100
    if i > 1        
        elimInd= cat(2,elimInd, elimInd+((length(Shuffle.mrall)/100)*(i-1)));
    end

end
newShuffle = Shuffle;
newShuffle.mrall(elimInd) = [];
newShuffle.firing(:,elimInd) = [];
end