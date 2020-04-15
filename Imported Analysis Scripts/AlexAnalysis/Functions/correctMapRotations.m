function returnP = correctMapRotations(paths)
    
    clc
    fprintf('\n')
    %%% Reliability constaint
    
    %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
    
    returnP = [];
    for mi = 1:length(upiece)
        fprintf(['\n\n\tMouse:  ' num2str(upiece{mi}) '\n'])
        
        isM = find(ismember(piece,upiece(mi)));
        s = load(paths{isM(1)});
        alignMap = s.alignment.alignmentMap;
        s.alignment.rotation = nan(size(alignMap));
        s.alignment.pval = repmat({[]},(size(alignMap)));
        redo = true;
        while redo
            for di = 1 %:1:length(isM);
                for dj = di+1:1:length(isM);

                    fprintf(['\n\t\tSessions:  ' num2str(spiece{isM(di)}) '\t' num2str(spiece{isM(dj)})])
                    if isempty(alignMap{di,dj})
                        if dj == length(isM)
                            redo = false;
                        end
                        continue
                    end
                    tic
                    s1 = load(paths{isM(di)});
                    s2 = load(paths{isM(dj)});

                    isInRoom1 = isInROI(s1.processed.p,s1.processed.roi.room);
                    isInRoom2 = isInROI(s2.processed.p,s2.processed.roi.room);

                    P1 = s1.processed.p;
                    P2 = s2.processed.p;
                    T1 = s1.processed.trace(alignMap{di,dj}(:,1),:);
                    T2 = s2.processed.trace(alignMap{di,dj}(:,2),:);

                    T = [T1 T2];

                    %%% Fit the rotation
                    doRot = [0:90:270];
                    bv = nan(size(doRot));
                    for r = doRot
                        RP2 = [cosd(r) -sind(r); sind(r) cosd(r)]*P2;
                        RP2 = bsxfun(@minus,RP2,nanmin(RP2,[],2));

                        [map samp val] = getMatchedSamplingMaps([P1 RP2],T,...
                            [true(1,length(P1)) false(1,length(P2))],...
                            [false(1,length(P1)) true(1,length(P2))]);
                        bv(r==doRot) = val;
                    end
                    [a b] = nanmax(bv);

                    if b~=1
                        fprintf('\n')
                        RP2 = [cosd(doRot(b)) -sind(doRot(b)); sind(doRot(b)) cosd(doRot(b))]*P2;
                        RP2 = bsxfun(@minus,RP2,nanmin(RP2,[],2));
                        s2.processed.p = RP2;
                        save(paths{isM(dj)},'-struct','s2','-v7.3');
                        returnP = [returnP; paths(isM(dj))];
                        break
                    elseif dj==length(isM)
                        redo = false;
                    end                
                end
            end
            
        end
    end
end












