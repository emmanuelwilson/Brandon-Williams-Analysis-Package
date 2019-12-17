%%Will take all of your MAT files and will create a single TIFF
%you will need the Multipage TIFF stack package 

function prepForCNMFE(folder)

name = dir(folder);
clipNum = nan(1,length(name));
for i = 1:length(name)
    if name(i).isdir == 0
        clipNum(i) = str2num(name(i).name(12:end-4));
    end
end
[a b] = sort(clipNum);
name = name(b);
clipNum = clipNum(b);

allFrames = [];


for i = 1:length(name)
    if ~isnan(clipNum(i))
        fprintf([num2str(clipNum(i)) '... ']);
        tmp = load([folder(1:end) '\' name(i).name]);
        tmp = tmp.Mr;
        %         tmp = uint8(tmp.Mr);
        %         tmp = double(tmp.Mr);
        
        if i == 1
            options.overwrite = true;
            options.append = false;
            options.big = true;
            saveastiff(tmp,[folder(1:end-19) '\' 'compiled.tif'],options);
        else
            options.big = true;
            options.overwrite = false;
            options.append = true;
            saveastiff(tmp,[folder(1:end-19) '\' 'compiled.tif'],options);
        end
    end
end
end