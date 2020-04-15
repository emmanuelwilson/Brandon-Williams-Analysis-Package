function ephysShift(paths)
    crossAllComp = [];
    crossRateComp = [];
    allRotCorrs = [];
    for p = paths'
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
%         gT = s.processed.trace(s.processed.splithalf.p<=0.05,:);
        gT = s.processed.trace;

        fprintf(['\n\t' p{1} ':  ' num2str(length(gT(:,1)))])

        [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        a = nanmin(indexSinceIn).*(1./30);
%         thresh = median(a(isInRoom));
        thresh = -1;
        
        indexSinceIn = nanmin(indexSinceIn).*(1./30) >= thresh;
        
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
        
        for door = 1:2
            comMaps(:,:,:,door) = mkTraceMaps(s.processed.p(:,isInRoom),...
                gT(:,isInRoom),isMostRecent(door,isInRoom) & half(1,isInRoom) & indexSinceIn(1,isInRoom), ...
                ones(1,2).*nanmax(doSize(1:2)));
            comMaps(:,:,:,door+2) = mkTraceMaps(s.processed.p(:,isInRoom),...
                gT(:,isInRoom),isMostRecent(door,isInRoom) & ~half(1,isInRoom) & indexSinceIn(1,isInRoom), ...
                ones(1,2).*nanmax(doSize(1:2)));
            
            mfr(:,door) = nanmean(gT(:,isMostRecent(door,:) & isInRoom & half),2);
            mfr(:,door+2) = nanmean(gT(:,isMostRecent(door,:) & isInRoom & ~half),2);

            
            totalMaps(:,:,:,door) = mkTraceMaps(s.processed.p(:,isInRoom),...
                gT(:,isInRoom),isMostRecent(door,isInRoom), ...
                ones(1,2).*nanmax(doSize(1:2)));
        end
        
%         for i = 1:4
%             
%         end
%         figure(1)
%         imagesc(pvxcorr(comMaps(:,:,20,3),comMaps(:,:,20,3),[10 10]))
%         colormap jet
%         axis equal
%         axis off
        
        
        
        doK = [3 3];
        
        for part = 0:floor(length(totalMaps(1,1,:,1))/prod(doK))
            
            figure(1)
            set(gcf,'position',[50 50 1350 900])
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(gT(:,1))
                    break
                end
                tmp = [totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,2) ...
                    nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,3)];
            
                subplot(doK(1),doK(2),k)
%                 
%                 step = 3;
%                 plot(s.processed.p(2,1:step:end),s.processed.p(1,1:step:end), ...
%                     'color',[0.0 0.0 0.0],'linewidth',1) %[0.8 0.8 0.8]
%                 hold on
%                 
%                 plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))&isMostRecent(1,:)),...
%                     s.processed.p(1,logical(gT(part.*prod(doK)+k,:))&isMostRecent(1,:)),...
%                     'color',[0.5 0.1 1],'linestyle','none','marker','o',...
%                     'markerfacecolor',[0.5 0.1 1],'markersize',2);
% 
%                 plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))&isMostRecent(2,:)),...
%                     s.processed.p(1,logical(gT(part.*prod(doK)+k,:))&isMostRecent(2,:)),...
%                     'color',[0.1 1 1],'linestyle','none','marker','o',...
%                     'markerfacecolor',[0.1 1 1],'markersize',2);

                imagesc(tmp)
                colormap jet
                alpha(double(~isnan(tmp)))
                set(gca,'ydir','normal')
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
        
%         crossRateComp = cat(3,crossRateComp,rateComp);
%         crossAllComp = cat(3,crossAllComp,allComp);
%         allRotCorrs = cat(4,allRotCorrs,rotCorrs);
    end
    crossAllComp = crossAllComp(1:2,3:4,:);
%     crossRateComp = crossRateComp(1:4,5:8,:);
    
%     clear v
%     mask = false(4,4);
%     mask(1:5:end) = true;
%     for i = 0:3
%         v(:,i+1) = help_getMaskedVals(crossRateComp,circshift(mask,[0 i]));
%     end
%     
%     figure(1)
%     set(gcf,'position',[50 50 250 250])
%     cumHist(v,[0:0.001:0.1])
%     axis square
%     
%     pval = nan(4,4);
%     for i = 1:4
%         for j = i+1:4
%             [h p ci tstat] = ttest(v(:,i),v(:,j));
%             pval(i,j) = p;
%         end
%     end
% %     imagesc(nanmean(crossRateComp(:,:,1:end),3))
%     
%     crossAllComp = cat(3,nanmean(crossAllComp(:,:,1:6),3), ... 
%         nanmean(crossAllComp(:,:,7:12),3), ...
%         nanmean(crossAllComp(:,:,13:18),3));
    
%     crossAllComp = cat(3,nanmean(crossAllComp(:,:,1:6),3),...
%         nanmean(crossAllComp(:,:,7:end),3));

    %%% Average within animal
    piece = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
    end
    upiece = unique(piece);
    val = [];
%     figure(3)
%     hold on
    for i = 1:length(upiece)
        v = nan(nansum(ismember(piece,upiece(i))),2);
        mask = false(2,2);
        mask(1:3:end) = true;
        v(:,1) = help_getMaskedVals(crossAllComp(:,:,ismember(piece,upiece(i))),mask);

        mask = false(2,2);
        mask(1,2) = true;
        mask(2,1) = true;
        v(:,2) = help_getMaskedVals(crossAllComp(:,:,ismember(piece,upiece(i))),mask);

%         plot(v(:,1)-v(:,2),'color','k');

        val = cat(3,val,nanmean(crossAllComp(:,:,ismember(piece,upiece(i))),3));
    end
    crossAllComp = val;

    %%%%%%%%%
    
    close all
    figure(1)
    set(gcf,'position',[50 50 425 300])
    imagesc(nanmean(crossAllComp(:,:,1:end),3))
    colormap(circshift([linspace(0,1,256)' ...
        [linspace(0,1,128) ones(1,128)]' ...
        linspace(0,1,256)'],[0 -2]))
    caxis([0.0 0.5])
    colorbar
    axis equal
    axis off
    
    
    
    allSim = crossAllComp;
    
    v = nan(length(allSim(1,1,:)),2);
    mask = false(2,2);
    mask(1:3:end) = true;
    v(:,1) = help_getMaskedVals(allSim,mask);
    
    mask = false(2,2);
    mask(1,2) = true;
    mask(2,1) = true;
    v(:,2) = help_getMaskedVals(allSim,mask);
    
    figure(2)
    set(gcf,'position',[50 450 250 200])
    mkGraph(v)
    
    
    
    root = 'Plots/Summary/TwoBig/';
    figure(1)
    saveFig(gcf,[root 'RDM'],'pdf');
    figure(2)
    saveFig(gcf,[root 'Similarity'],'pdf');
    
    outP = ['Stats_TwoBig.txt'];
    fid = fopen(outP,'w');
    fprintf(fid,'\t\t\tSIMILARITY\n');
    for i = 1:2
        for j = i+1:2
            [h p ci tstat] = ttest(v(:,i),v(:,j));
            fprintf(fid,['\n' num2str(i) ' to ' num2str(j) ':  ']);
            fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
        end
    end
    
    fclose all
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end









































