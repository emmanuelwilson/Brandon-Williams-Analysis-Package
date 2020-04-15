
function multiSessionReliability(paths)
    
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
                
                [isIn mrd1] = isInROI(s1.processed.p,(s1.processed.roi.door));
                [isIn mrd2] = isInROI(s2.processed.p,(s2.processed.roi.door));
                
                P1 = s1.processed.p(:,isInRoom1);
                P2 = s2.processed.p(:,isInRoom2);
                T1 = s1.processed.trace(alignMap{di,dj}(:,1),isInRoom1);
                T2 = s2.processed.trace(alignMap{di,dj}(:,2),isInRoom2);
                P = [P1 P2];
                T = [T1 T2];
                
                allMasks = repmat({[]},[1 length(s.processed.roi.door(1,:)).*2]);
                for doorA = 1:length(s.processed.roi.door(1,:))
                    allMasks{doorA} = [mrd1(doorA,isInRoom1) false(1,length(mrd2(doorA,isInRoom2)))];
                    allMasks{doorA+length(s.processed.roi.door(1,:))} = ...
                        [false(1,length(mrd1(doorA,isInRoom1))) mrd2(doorA,isInRoom2)];
                end
                
                [map samp val ival] = getMatchedMapsNMasks(P,T,allMasks);
                tmp = ival(1:length(s.processed.roi.door(1,:)),...
                    length(s.processed.roi.door(1,:))+1:end,:);
                tmp(isnan(tmp)) = -inf;
                tmp = sort(reshape(tmp,numel(tmp(:,:,1)),[]),'descend');
                actual = nanmean(tmp(1:length(s.processed.roi.door(1,:)),:))';

                nsims = 100;
                minShift = 900;
                null = nan(length(actual),nsims);
                fprintf(['\n\t\t\t( ' num2str(di) ', ' num2str(dj) '; rotation:  ' ...
                    num2str(s.alignment.rotation(di,dj)) ')  Computing null...'])
                parfor sim = 1:nsims
                    sp = circshift(P,[0 minShift+randi(length(P)-minShift.*2)]);
                    [map samp val ival] = getMatchedMapsNMasks(sp,T,allMasks);
                    
                    
                    tmp = ival(1:length(s.processed.roi.door(1,:)),...
                        length(s.processed.roi.door(1,:))+1:end,:);
                    tmp(isnan(tmp)) = -inf;
                    tmp = sort(reshape(tmp,numel(tmp(:,:,1)),[]),'descend');
                    null(:,sim) = nanmean(tmp(1:length(s.processed.roi.door(1,:)),:))'
                end
                pvals = 1-nanmean(bsxfun(@gt,actual',null(:)))';
                s.alignment.pval{di,dj} = pvals;
                
                close all
                set(gcf,'position',[600 50 250 250])
                cumHist({pvals},[0:0.01:1]);
                drawnow
                
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















