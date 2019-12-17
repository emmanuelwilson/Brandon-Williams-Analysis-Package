function out = PokeFind(sep)
fnames = fieldnames(sep);
if any(strcmp(fnames, 'CorrectTrial'))
c = find(~cellfun(@isempty,sep.CorrectTrial));
else
    c = [];
end
if any(strcmp(fnames, 'IncorrectTrial'))
ic = find(~cellfun(@isempty,sep.IncorrectTrial));
else
    ic = [];
end
if any(strcmp(fnames, 'CorrectCorrectionTrial'))
cc = find(~cellfun(@isempty,sep.CorrectCorrectionTrial));
else
    cc = [];
end
if any(strcmp(fnames, 'IncorrectCorrectionTrial'))
icc = find(~cellfun(@isempty,sep.IncorrectCorrectionTrial));
else
    icc = [];
end
out = unique(cat(2,c,ic,cc,icc));
end