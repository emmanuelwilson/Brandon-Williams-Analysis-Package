function catMiniscopeFiles(allPaths,outFolder)
    if nargin<2
        outFolder = [inFolder '/Merged'];
    end

    checkP([outFolder '/blah']);
    
    doExclude = [];
    newTs = [];
    currentTs = 0;
    currentBV = 0;
    currentMsV = 0;
    currentBVFrameCount = 0;
    currentMsVFrameCount = 0;
    partitions = [];
    for i = 1:length(allPaths)       
        path = allPaths{i};
        slashInds = find(ismember(path,'/'));
        bv = dir([path(1:slashInds(end)) 'behavCam*.avi']);
        bv = cat(2,{bv(:).name});
        
        totalBVFrames = 0;
        for bvi = 1:length(bv)
            v = VideoReader([path(1:slashInds(end)) bv{bvi}]);
            totalBVFrames = totalBVFrames + v.NumberOfFrames;
        end
        
        fid = fopen(allPaths{i});
        dataArray = textscan(fid, '%f%f%f%f%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
        fclose(fid);
        isGood = false(length(dataArray{1}(:,1)),1);
        isGood(find(abs(dataArray{3}(1:15)) > 500)) = true;
        doExclude = [doExclude; isGood];
        og = cell2mat(dataArray(1:4));
        
        if nansum(og(:,1)==1) == totalBVFrames
            behavCam = 1;
            msCam = 0;
        else
            behavCam = 0;
            msCam = 1;
        end
        
        tmp = og;
        tmp(:,3) = tmp(:,3) + currentTs;
        tmp(tmp(:,1)==behavCam,2) = tmp(tmp(:,1)==behavCam,2) + currentBVFrameCount;
        tmp(tmp(:,1)==msCam,2) = tmp(tmp(:,1)==msCam,2) + currentMsVFrameCount;
        
        newTs = [newTs; tmp(:,1:2) tmp(:,3)+currentTs tmp(:,4)];
        
        currentTs = currentTs + tmp(end,3);
        partitions = [partitions currentTs];
        
        blah = tmp(tmp(:,1)==behavCam,2);
        currentBVFrameCount = blah(end);
        blah = tmp(tmp(:,1)==msCam,2);
        currentMsVFrameCount = blah(end);
        
        clipNum = nan(1,length(bv));
        for i = 1:length(bv)
            clipNum(i) = str2num(bv{i}(9:end-4));
        end
        [a b] = sort(clipNum);
        bv = bv(b);
        clipNum = clipNum(b);
        
        for bvi = 1:length(bv)
            copyfile([path(1:slashInds(end)) bv{bvi}],[outFolder '/behavCam' num2str(bvi+currentBV) '.avi']);
        end
        
        currentBV = currentBV+length(bv);
        
        msv = dir([path(1:slashInds(end)) 'msCam*.avi']);
        msv = cat(2,{msv(:).name});
        
        clipNum = nan(1,length(msv));
        for i = 1:length(msv)
            clipNum(i) = str2num(msv{i}(6:end-4));
        end
        [a b] = sort(clipNum);
        msv = msv(b);
        clipNum = clipNum(b);
        
        for msvi = 1:length(msv)
            copyfile([path(1:slashInds(end)) msv{msvi}],[outFolder '/msCam' num2str(msvi+currentMsV) '.avi']);
        end
        currentMsV = currentMsV+length(msv);
    end
    tsTitle = ['camNum	frameNum	sysClock	buffer\n'];
    fid = fopen([outFolder '/timestamp.dat'],'w');
    fprintf(fid,tsTitle);
    for i = 1:length(newTs)
        fprintf(fid,'%i\t%i\t%i\t%i\n',newTs(i,:));
    end
    fclose(fid);
    
    fid = fopen([outFolder '/exclude.dat'],'w');
    fprintf(fid,tsTitle);
    for i = 1:length(doExclude)
        fprintf(fid,'%i\n',doExclude(i));
    end
    fclose(fid);
    
    fid = fopen([outFolder '/partitions.dat'],'w');
    for i = 1:length(partitions)
        fprintf(fid,'%i\n',partitions(i));
    end
    fclose(fid);
end


















