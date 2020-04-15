function prepForCNMFE(folder)
    
    name = dir(folder);
    clipNum = nan(1,length(name));
    for i = 1:length(name)
        clipNum(i) = str2num(name(i).name(12:end-4));
    end
    [a b] = sort(clipNum);
    name = name(b);
    clipNum = clipNum(b);
    
    allFrames = [];
    
    
    for i = 1:length(name)
        fprintf([num2str(clipNum(i)) '... ']);
        tmp = load([folder(1:end-5) name(i).name]);
        tmp = tmp.Mr;
%         tmp = uint8(tmp.Mr);
%         tmp = double(tmp.Mr);
        
        if i == 1
            options.overwrite = true;
            options.append = false;
            options.big = true;
            saveastiff(tmp,[folder(1:end-19) 'compiled.tif'],options);
        else
            options.big = true;
            options.overwrite = false;
            options.append = true;
            saveastiff(tmp,[folder(1:end-19) 'compiled.tif'],options);
        end
    end
end