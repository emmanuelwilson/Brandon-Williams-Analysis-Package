%% Will Eliminate any cells that do not have a peak within 60 frames of reward
function [out] = FindRewardCells(ShuffledCrit)
CellNum = ShuffledCrit.BcorrectNumberPassed;
CellInd = ShuffledCrit.BcorrectCells;
Traces = permute(ShuffledCrit.BcorrectPassed,[3 2 1]);
ITI = 450;
graceP = 60;

for i = CellNum: -1 : 1 
    [maxvals,mtimes] = max(Traces(:,:,i),[],2);
    if length(mtimes) == length(find(any(mtimes < (ITI-graceP)))) || length(mtimes) == length(find(any(mtimes > (ITI+graceP))))
        CellNum = CellNum-1;
        CellInd(i) = [];
        Traces(:,:,i) = [];
    end
end
meanFluo = nanmean(Traces,1);
out.PopulationFluo = permute(meanFluo,[3 2 1]);
out.NumberOfCells = CellNum;
out.CellInd = CellInd;
out.Traces = Traces; 
end