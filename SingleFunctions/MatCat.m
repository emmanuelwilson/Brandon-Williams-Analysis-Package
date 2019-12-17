%%This function will take all the alligned MAT files and concatenate them
function MatCat(folder)


folder = [folder '/Normcorred*.mat'];
%%% Sort Clips in chronological order
name = dir(folder);
clipNum = nan(1,length(name));
for i = 1:length(name)
    clipNum(i) = str2num(name(i).name(12:end-4));
end
[a b] = sort(clipNum);
name = name(b);
clipNum = clipNum(b);
Mrtot = [];

for i = 1 : length(clipNum)
    if i ==1
        load([folder(1:end-16) '\' name(i).name])
        Mrtot = Mr;
    else
        load([folder(1:end-16) '\' name(i).name])
        Mrtot = cat(3,Mrtot,Mr);
    end
end
Y = double(Mrtot);
Ysiz = size(Y);
save('AlignedMAT.mat','Y','Ysiz', '-v7.3')
end