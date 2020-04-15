
function multiSessionReliability(paths)
    
    clc
    fprintf('\n')
    %%% Reliability constaint
    
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
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
    
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}) '\n'])
        
        isM = find(ismember(piece,upiece(mi)));
        s = load(paths{isM(1)});
        alignMap = s.alignment.alignmentMap;
        s.alignment.rotation = nan(size(alignMap));
        s.alignment.pval = repmat({[]},(size(alignMap)));
        for di = 1:1:length(isM);
            for dj = di+1:1:length(isM);
                
                fprintf(['\n\t\tSessions:  ' num2str(spiece{isM(di)}) '\t' num2str(spiece{isM(dj)})])
                if isempty(alignMap{di,dj})
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
                doRot = [0];
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
                s.alignment.rotation(di,dj) = doRot(b);
                
                RP2 = [cosd(doRot(b)) -sind(doRot(b)); sind(doRot(b)) cosd(doRot(b))]*P2;
                RP2 = bsxfun(@minus,RP2,nanmin(RP2,[],2));
                P = [P1 RP2]; %%% Use optimal rotation
                
                [map samp val ival] = getMatchedSamplingMaps(P,T,...
                    [true(1,length(P1)) false(1,length(P2))],...
                    [false(1,length(P1)) true(1,length(P2))]);
                
                actual = ival;

                nsims = 100;
                minShift = 900;
                null = nan(length(actual),nsims);
                fprintf(['\n\t\t\t( ' num2str(di) ', ' num2str(dj) '; rotation:  ' ...
                    num2str(s.alignment.rotation(di,dj)) ')  Computing null...'])
                parfor sim = 1:nsims
                    sp = circshift(P,[0 minShift+randi(length(P)-minShift.*2)]);
                    [map samp val ival] = getMatchedSamplingMaps(sp,T,...
                        [true(1,length(P1)) false(1,length(P2))],...
                        [false(1,length(P1)) true(1,length(P2))]);
                    null(:,sim) = ival;
                end
                pvals = 1-nanmean(bsxfun(@gt,actual',null(:)))';
                s.alignment.pval{di,dj} = pvals;
                
                durat = toc;
                fprintf(['\n\t\t\tCount:  ' num2str(nansum(pvals <= 0.05)) ...
                    '\tProportion:  ' num2str(nanmean(pvals <= 0.05)) '\tDuration:  ' num2str(durat) '\n'])

            end
        end
        save(paths{isM(1)},'-struct','s','-v7.3');
    end
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end















