%When you forget to eliminate duplicate cells before shuffling

function [out] = ElimdupfromShuffle(Shuf,elim)
indOut = [];
for i = 1: length(elim)   
    indOut(i,:) = [elim(i):length(Shuf.mrall)/100:length(Shuf.mrall)];
end

indOut = sort(reshape(indOut, [1, length(indOut(:,1))*length(indOut(1,:))]),'descend');
out = Shuf;
for i = 1 : length(indOut)
    out.mrall(indOut(i)) = [];
% out.firing(:,elim) = [];
end
out.percentil99th = prctile(out.mrall,99);
end