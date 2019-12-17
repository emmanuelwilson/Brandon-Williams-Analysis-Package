function [corrmat,pairs] = duplicateFinder(ms)

thresh = 0.6;

corrmat = NaN(length(ms.FiltTraces(1,:)),length(ms.FiltTraces(1,:)));

for i = 1 : length(corrmat(:,1))
    for j = 1 : length(corrmat(:,1))
        corrmat(i,j) = corr2(ms.FiltTraces(:,i),ms.FiltTraces(:,j));
    end
end

figure
imagesc(corrmat)
colorbar

[row, col] = find(corrmat>thresh & corrmat<1);
pairs = cat(2,col,row);
dup = [];
count = 1;
for i = 1: length(pairs(:,1))
    if pairs(i,1)>pairs(i,2) || find(pairs(:,2) == pairs(i,1),1,'first')<i
        dup(count) = i;
        count = count +1;
    end
end
pairs(dup,:) = [];

end