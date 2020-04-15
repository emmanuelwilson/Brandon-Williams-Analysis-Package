function mapsXDoors(paths)
    crossAllComp = [];
    for p = paths'
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
        gT = s.processed.trace(s.processed.splithalf.p<=0.05,:);
%          gT = s.processed.trace(s.processed.splithalf.within.p<=0.05,:);
%         gT = s.processed.trace;

        fprintf(['\n\t' p{1} ':  ' num2str(length(gT(:,1)))])

        [isIn isMostRecent] = isInROI(s.processed.p,(s.processed.roi.door));
%         isMostRecent = fliplr(isMostRecent);
%         isMostRecent = nthMostRecent(isIn,2);
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        
        half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
        
        allComp = nan(length(s.processed.roi.door(1,:)).^2);
        for doorA = 1:length(s.processed.roi.door(1,:))
            for doorB = 1:length(s.processed.roi.door(1,:))
                for halfA = 0:1
                    for halfB = 0:1
                        if (halfA).*length(s.processed.roi.door(1,:))+doorA > ...
                                (halfB).*length(s.processed.roi.door(1,:))+doorB || halfA == halfB
                            continue
                        end
                        [map samp val] = getMatchedSamplingMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),...
                            [isMostRecent(doorA,isInRoom) & halfA==half(1,isInRoom)],...
                            [isMostRecent(doorB,isInRoom) & halfB==half(1,isInRoom)]);
                        allComp((halfA).*length(s.processed.roi.door(1,:))+doorA,...
                            (halfB).*length(s.processed.roi.door(1,:))+doorB) = val;
                    end
                end
            end
        end
        
        doSize = zeros(1,3);
        for room = 1:1
            m1 = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]));
            doSize = nanmax([doSize; size(m1)]);
        end
        
        comMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 8]);
        totalMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 4]);
        mfr = nan([doSize(3) 8]);
        
        half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
        
        for door = 1:length(s.processed.roi.door(1,:))
            comMaps(:,:,:,door) = mkTraceMaps(s.processed.p(:,isInRoom),...
                gT(:,isInRoom),isMostRecent(door,isInRoom) & half(1,isInRoom), ...
                ones(1,2).*nanmax(doSize(1:2)));
            comMaps(:,:,:,door+length(s.processed.roi.door(1,:))) = mkTraceMaps(s.processed.p(:,isInRoom),...
                gT(:,isInRoom),isMostRecent(door,isInRoom) & ~half(1,isInRoom), ...
                ones(1,2).*nanmax(doSize(1:2)));
            
            mfr(:,door) = nanmean(gT(:,isMostRecent(door,:) & isInRoom & half),2);
            mfr(:,door+length(s.processed.roi.door(1,:))) = nanmean(gT(:,isMostRecent(door,:) & isInRoom & ~half),2);

            
            totalMaps(:,:,:,door) = mkTraceMaps(s.processed.p(:,isInRoom),...
                gT(:,isInRoom),isMostRecent(door,isInRoom), ...
                ones(1,2).*nanmax(doSize(1:2)));
        end
       
        
        doK = [8 8];
        
        for part = 0:floor(length(totalMaps(1,1,:,1))/prod(doK))
            
            figure(1)
            set(gcf,'position',[50 50 900 900])
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(gT(:,1))
                    break
                end
                tmp = [[totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,2)] ; ...
                    nan(1,size(totalMaps,2).*2 + 1); ...
                    [totalMaps(:,:,part.*prod(doK)+k,3) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,4)]];
            
                subplot(doK(1),doK(2),k)
                imagesc(tmp)
                colormap jet
                alpha(double(~isnan(tmp)))
                axis equal
                axis off    
            end
            
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/DifferentiatedMaps/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'tiff')
            saveFig(gcf,outP,'pdf')
            close all
            drawnow
        end
        
        crossAllComp = cat(3,crossAllComp,allComp);
    end
    crossAllComp = crossAllComp(1:length(s.processed.roi.door(1,:)),...
        length(s.processed.roi.door(1,:))+1:length(s.processed.roi.door(1,:)).*2,:);

    %% Average within animal
    piece = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
    end
    upiece = unique(piece);
    val = [];
    for i = 1:length(upiece)
        val = cat(3,val,nanmean(crossAllComp(:,:,ismember(piece,upiece(i))),3));
    end
    crossAllComp = val;
    
    close all
    figure(1)
    set(gcf,'position',[50 50 425 300])
    imagesc(nanmean(crossAllComp(:,:,1:end),3))
    colormap(circshift([linspace(0,1,256)' ...
        [linspace(0,1,128) ones(1,128)]' ...
        linspace(0,1,256)'],[0 -2]))
    caxis([0.0 0.4])
    colorbar
    axis equal
    axis off
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end









































