function plotMultiDayMaps(paths)
    
    warning off all
    
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
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        alignMap = s.alignment.alignmentMap;
        am = repmat({[]},[1 length(sessions)]);
        isPC = repmat({[]},[1 length(sessions)]);
        for si = 1:length(sessions)
            s = load(sessions{si});
            slashInds = find(ismember(sessions{si},'/'));
            fprintf(['\t\tPreloading:  ' sessions{si}(slashInds(end)+1:end-4) '\n'])
            
%             [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
%             [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
%             allMasks = [{isMostRecent(1,isInRoom)}  {isMostRecent(2,isInRoom)}];
            gT = s.processed.trace;
%             mapA = mkTraceMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks{1},[13 13]);
%             mapB = mkTraceMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks{2},[13 13]);
            map2 = mkTraceMaps(s.processed.p,gT,[],[17 17]);
            am{si} = map2;
            isPC{si} = s.processed.splithalf.wholemap.p <= 0.025;
        end  
        
        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                tmp1 = am{si};
                tmp2 = am{sj};
                
                if isempty(alignMap{si,sj})
                    continue
                end
                
                
                isGood = isPC{si}(alignMap{si,sj}(:,1)) | isPC{sj}(alignMap{si,sj}(:,2));
                tmp1 = tmp1(:,:,alignMap{si,sj}(isGood,1));
                tmp2 = tmp2(:,:,alignMap{si,sj}(isGood,2));
                vals = nan(1,4);
                for rot = 0:3
                    rtmp2 = imrotate(tmp2,rot.*90);
                    goodPixels = ~isnan(rtmp2(:,:,1))&~isnan(tmp1(:,:,1));
                    vals(rot+1) = corr(tmp1(repmat(goodPixels,[1 1 length(tmp1(1,1,:))])), ...
                        rtmp2(repmat(goodPixels,[1 1 length(tmp1(1,1,:))])));
                end
                [a bestRot] = nanmax(vals);
                tmp2 = imrotate(tmp2,(bestRot-1).*90);
                
                toPlot = [tmp1 nan(length(tmp1(:,1,1)),4,nansum(isGood)) tmp2];
                
                doK = [8 4];

                for part = 0:floor(length(toPlot(1,1,:,1))/prod(doK))

                    figure(1)
                    set(gcf,'position',[50 50 900 1350])
                    for k = 1:prod(doK)
                        if part.*prod(doK)+k > length(toPlot(1,1,:))
                            break
                        end
                        
                        subplot(doK(1),doK(2),k)
                        imagesc(toPlot(:,:,part.*prod(doK)+k))
                        colormap jet
%                         caxis([0 nanmax(nanmax(toPlot(:,:,part.*prod(doK)+k)))])
                        alpha(double(~isnan(toPlot(:,:,part.*prod(doK)+k))))
                        axis equal
                        axis off    
                    end
                    slashInds1 = find(ismember(sessions{si},'/'));
                    slashInds2 = find(ismember(sessions{sj},'/'));
                    outP = ['Plots/PairwiseAlignedCellMaps/' upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end) '/' ...
                        [sessions{si}(slashInds1(end)+1:end-4) '_vs_' sessions{sj}(slashInds2(end)+1:end-4)] ...
                        '_Part_' num2str(part)];
                    saveFig(gcf,outP,[{'tiff'} {'pdf'}])
                    close all
                    drawnow
                end
            end
        end
    end
end