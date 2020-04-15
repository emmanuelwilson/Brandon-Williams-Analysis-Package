function s = alignTraceData(doPath,s)
    fid = fopen([doPath '/timestamp.dat']);
    dataArray = textscan(fid, '%f%f%f%f%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
    fclose(fid);    
    
    camera = dataArray{1};
    frame = dataArray{2};
    sysclock = dataArray{3};
        
    if nanmax(frame(camera==0)) == length(s.pos.mp(1,:,1)) & ...
            nanmax(frame(camera==1)) == length(s.calcium.trace(:,1))
        posCam = 0;
        traceCam = 1;
    elseif nanmax(frame(camera==1)) == length(s.pos.mp(1,:,1)) & ...
            nanmax(frame(camera==0)) == length(s.calcium.trace(:,1))
        posCam = 1;
        traceCam = 0;
    else
        fprintf('Frame time alignment error, cannot assign cameras')
    end
    
    if exist([doPath '/exclude.dat'],'file')==2
        fid = fopen([doPath '/exclude.dat']);
        exclude = textscan(fid, '%u', 'Delimiter', '\n', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
        fclose(fid);    
        camera(logical(exclude{1})) = [];
        frame(logical(exclude{1})) = [];
        sysclock(logical(exclude{1})) = [];
    else
        camera(find(abs(sysclock(1:15)) > 500)) = [];
        frame(find(abs(sysclock(1:15)) > 500)) = [];
        sysclock(find(abs(sysclock(1:15)) > 500)) = [];
    end
    
%     if exist([doPath '/partitions.dat'],'file')==2
%         fid = fopen([doPath '/partitions.dat']);
%         partitions = textscan(fid, '%u', 'Delimiter', '\n', 'EmptyValue' ,NaN,'HeaderLines' ,0, 'ReturnOnError', false);
%         fclose(fid);    
%     end

    
%     oldP = s.pos.p;
%     oldP = bsxfun(@minus,oldP,nanmin(oldP')');
%     oldP = bsxfun(@times,oldP,nanmax(s.environment.size./nanmax(oldP')));
%     
%     
%     unP = s.pos.uninterp;
%     unP = bsxfun(@minus,unP,nanmin(unP')');
%     unP = bsxfun(@times,unP,nanmax(s.environment.size./nanmax(unP')));
%     [unP gaps] = interpNaNs(unP');
%     unP = unP';
%     unP(:,gaps>=15) = nan;
    
    posFrames = [frame(camera==posCam) sysclock(camera==posCam)];
    traceFrames = [frame(camera==traceCam) sysclock(camera==traceCam)];

    traceFrames(traceFrames(:,2) >= posFrames(end,2),:) = [];
    traceFrames(traceFrames(:,2) <= posFrames(1,2),:) = [];
%     
%     linX = linterp(posFrames(:,2),oldP(1,posFrames(:,1)),traceFrames(:,2));
%     linY = linterp(posFrames(:,2),oldP(2,posFrames(:,1)),traceFrames(:,2));
    
    for i = 1:length(s.pos.mp(1,1,:))
        unX(:,:,i) = linterp(posFrames(:,2),s.pos.mp(1,posFrames(:,1),i),traceFrames(:,2));
        unY(:,:,i) = linterp(posFrames(:,2),s.pos.mp(2,posFrames(:,1),i),traceFrames(:,2));
    end
    
%     s.processed.p = [linX'; linY'];
    s.processed.mp = [permute(unX,[2 1 3]); permute(unY,[2 1 3])];
    s.processed.trace = s.calcium.trace(traceFrames(:,1),:)';
    s.processed.validTraceFrames = traceFrames;
    s.processed.posFrames = posFrames;
%     s.properties.partitions = [find(diff(posFrames(:,2))>300) length(posFrames(:,2))];
end