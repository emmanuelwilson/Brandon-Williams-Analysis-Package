function mapsXDoors(paths)
    crossAllComp = [];
    allRotCorrs = [];
    for p = paths'
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
%         gT = s.processed.trace(s.processed.splithalf.p<=0.05,:);
        gT = s.processed.trace;
        
        [isIn isMostRecent] = isInROI(s.processed.p,s.processed.roi.door);
        
        com = isMostRecent(3,:);
        uni = isMostRecent(1,:)|isMostRecent(2,:);
        
        isInRoom = isInROI(s.processed.p,s.processed.roi.room);
        
        
        doSize = zeros(1,3);
        for room = 1:2
            m1 = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),com(isInRoom(room,:)));
            m2 = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),uni(isInRoom(room,:)));
            doSize = nanmax([doSize; size(m1); size(m2)]);
        end
        
        comMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 4]);
        uniMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 4]);
        
        totalMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 4]);
        
        half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
        for room = 1:2
            comMaps(:,:,:,room) = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),com(isInRoom(room,:)) & half(isInRoom(room,:)),ones(1,2).*nanmax(doSize(1:2)));
            uniMaps(:,:,:,room) = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),uni(isInRoom(room,:)) & half(isInRoom(room,:)),ones(1,2).*nanmax(doSize(1:2)));
            comMaps(:,:,:,room+2) = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),com(isInRoom(room,:)) & ~half(isInRoom(room,:)),ones(1,2).*nanmax(doSize(1:2)));
            uniMaps(:,:,:,room+2) = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),uni(isInRoom(room,:)) & ~half(isInRoom(room,:)),ones(1,2).*nanmax(doSize(1:2)));
            
            totalMaps(:,:,:,room+2) = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),com(isInRoom(room,:)),ones(1,2).*nanmax(doSize(1:2)));
            totalMaps(:,:,:,room) = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),uni(isInRoom(room,:)),ones(1,2).*nanmax(doSize(1:2)));
        end
        
        comMaps = cat(4,comMaps(:,:,:,1:2),uniMaps(:,:,:,1:2), ...
            comMaps(:,:,:,3:4),uniMaps(:,:,:,3:4));
        allComp = nan(8,8);
        rotCorrs = nan(4,2);
        for ri = 1:8
            for rj = 1:8
                m1 = comMaps(:,:,:,ri);
                m2 = comMaps(:,:,:,rj);
                if isempty(m1) || isempty(m2)
                    continue
                end
                
                isGood = ~(isnan(m1)|isnan(m2));
                if ~any(isGood(:))
                    continue
                end
                
                allComp(ri,rj) = corr(m1(isGood),m2(isGood));
            end
        end
        
%         doK = [8 8];
%         
%         for part = 0:floor(length(totalMaps(1,1,:,1))/prod(doK))
%             
%             figure(1)
%             set(gcf,'position',[50 50 900 900])
%             for k = 1:prod(doK)
%                 if part.*prod(doK)+k > length(gT(:,1))
%                     break
%                 end
%                 tmp = [[totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,3)] ; ...
%                     nan(1,size(totalMaps,2).*2 + 1); ...
%                     [totalMaps(:,:,part.*prod(doK)+k,2) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,4)]];
%             
%                 subplot(doK(1),doK(2),k)
%                 imagesc(tmp)
%                 colormap jet
%                 alpha(double(~isnan(tmp)))
%                 axis equal
%                 axis off    
%             end
%             
%             slashInds = find(ismember(p{1},'/'));
%             outP = ['Plots/DifferentiatedMaps/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
%             saveFig(gcf,outP,'tiff')
%             saveFig(gcf,outP,'pdf')
%             close all
%             drawnow
%         end

        
%         imagesc(rotCorrs)
        
        crossAllComp = cat(3,crossAllComp,allComp(1:4,5:8));
%         allRotCorrs = cat(3,allRotCorrs,rotCorrs);
    end
    
    close all
    figure(1)
    set(gcf,'position',[50 50 425 300])
    imagesc(nanmean(crossAllComp(:,:,7:end),3))
    colormap(circshift([linspace(0,1,256)' ...
        [linspace(0,1,128) ones(1,128)]' ...
        linspace(0,1,256)'],[0 -2]))
    caxis([-0.01 0.5])
    colorbar
    axis equal
    axis off
    
    allSim = crossAllComp;
    
    v = nan(length(allSim(1,1,:)),6);
    mask = false(4,4);
    mask(11:5:end) = true;
    v(:,1) = help_getMaskedVals(allSim,mask);
    
    mask = false(4,4);
    mask(1:5:10) = true;
    v(:,2) = help_getMaskedVals(allSim,mask);
    
    mask = false(4,4);
    mask(4,3) = true;
    mask(3,4) = true;
    v(:,3) = help_getMaskedVals(allSim,mask);
    
    mask = false(4,4);
    mask(2,1) = true;
    mask(1,2) = true;
    v(:,4) = help_getMaskedVals(allSim,mask);
    
    
    mask = false(4,4);
    mask(1:5:10) = true;
    v(:,6) = help_getMaskedVals(allSim,mask);
    mask(11:5:end) = true;
    v(:,5) = help_getMaskedVals(allSim,mask);
    mask = false(4,4);
    mask(1,3) = true;
    mask(3,1) = true;
    mask(2,4) = true;
    mask(4,2) = true;
    v(:,7) = help_getMaskedVals(allSim,mask);
    
    
    figure(2)
    set(gcf,'position',[50 450 800 200])
    toPlot = v';
    toPlot = num2cell(toPlot);
    subplot(1,6,1:2)
    mkGraph(toPlot([1 2],:)')
%     set(gca,'ylim',[0 0.7])
    subplot(1,6,3:4)
    mkGraph(toPlot([3 4],:)')
%     set(gca,'ylim',[-0.15 0.15])
    subplot(1,6,5)
    mkGraph(v(7:end,[5:7]))
    subplot(1,6,6)
%     mkGraph(v(7:end,[5:7])')
    mkGraph([nanmean(v(7:end,[5:6]),2)-v(7:end,[7])]')
%     set(gca,'ylim',[0 0.7])
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end












































