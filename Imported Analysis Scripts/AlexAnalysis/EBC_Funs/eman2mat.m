function eman2mat(inFolder)
    
    
    warning off all

    p = getFilePaths('MatlabData_Emmanuel','.mat');
    
    pieces = [];
    for i = 1:length(p)
        pieces = [pieces; {p{i}(1:find(ismember(p{i},'/'),1,'last')-1)}];
    end
    upiece = unique(pieces);
    failed = [];
    for i = 1:length(upiece)
        s = struct;
        slashInds = find(ismember(upiece{i},'/'));
        s.properties.mouse = upiece{i}(slashInds(1)+1:slashInds(2)-1);
        s.properties.session = upiece{i}(slashInds(2)+1:end);
        
        fprintf(['\n\t' s.properties.mouse '\t' s.properties.session]);
        
        outP = ['MatlabData/' s.properties.mouse '/' s.properties.session];
        if exist([outP '.mat'])==2
            continue
        end
        
        try

            p = getFilePaths(upiece{i},'msDeconvolved.mat');
            tmp = load(p{1});
            s.calcium = tmp.ms;
            s.processed.trace = s.calcium.FiltTraces;
            s.properties.date = [tmp.ms.dirName(end-5:end-4) '-' ...
                tmp.ms.dirName(end-3:end-2) '-' tmp.ms.dirName(end-1:end-0)];
            
            p = getFilePaths(upiece{i},'HeadTrackingData.mat');
            tmp = load(p{1});
            
            s.pos.p = tmp.SINKdata(:,1:2)';
            s.pos.hd = tmp.HDdeg';
            
            
            p = getFilePaths(upiece{i},'frameMap.mat');
            tmp = load(p{1});
            s.frameMap = tmp.frameMap;
            
            outP = ['MatlabData/' s.properties.mouse '/' s.properties.session];
            checkP(outP);
            save(outP,'-struct','s','-v7.3');
            
        catch
            failed = [failed; upiece(i)];
        end
    end    
end