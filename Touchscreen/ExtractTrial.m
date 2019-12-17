function [out,name] = ExtractTrial(rasts,num,choice)
if num == 1
    out = rasts.CorrectTrial{choice};    
elseif num == 2
    out = rasts.IncorrectTrial{choice};
elseif num ==3
    out = rasts.CorrectCorrectionTrial{choice};
elseif num == 4
    out = rasts.IncorrectCorrectionTrial{choice};
end
end