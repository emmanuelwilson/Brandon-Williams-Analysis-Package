function [dcoef, ndcoef] = DelayActivityCoef(trace,dstart,dFrames)
nondelay = trace;
delay = nan((dFrames+1)*length(dstart),1);
for i = length(dstart): -1: 1
    nondelay(dstart(i):dstart(i)+dFrames) = [];
    delay((i*(dFrames+1))-dFrames:(i*(dFrames+1))) = trace(dstart(i):dstart(i)+dFrames);
end

% d = Binarize(delay);
dcoef = length(find(delay))/length(delay);
% nd = Binarize(nondelay);
ndcoef = length(find(nondelay))/length(nondelay);
end