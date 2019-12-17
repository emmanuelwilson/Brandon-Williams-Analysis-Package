%%This Script will convert all of your "msCam.AVI" videos in your directory
%%into one large tiff file.

folder = [pwd '/Normcorred*.mat'];

name = dir(folder);
clipNum = nan(1,length(name));
for i = 1:length(name)
    clipNum(i) = str2num(name(i).name(12:end-4));
end
[a, b] = sort(clipNum);
name = name(b);
% name.name = [folder, '\', name(:).name;];
for i =1 : length(clipNum)
    load(strcat(folder(1:end-16),'\', name(i).name));
    Mr = uint16(Mr);
    frames = length(Mr(1,1,:));
    for x = 1 : frames
        if i == 1 && x ==1
            imwrite(Mr(:,:,x),'Aligned.tif');
        else
            imwrite(Mr(:,:,x),'Aligned.tif','WriteMode','append');
        end
    end
end
