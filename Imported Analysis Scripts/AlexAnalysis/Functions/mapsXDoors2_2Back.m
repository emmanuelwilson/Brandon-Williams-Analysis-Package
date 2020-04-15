function mapsXDoors2(paths)
    crossAllComp = [];
    crossRateComp = [];
    allRotCorrs = [];
    for p = paths'
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
        gT = s.processed.trace(s.processed.splithalf.p<=0.05,:);
%         gT = s.processed.trace;

        fprintf(['\n\t' p{1} ':  ' num2str(length(gT(:,1)))])

        [isIn isMostRecent] = isInROI(s.processed.p,(s.processed.roi.door));
%         isMostRecent = fliplr(isMostRecent);
        isSecondMostRecent = nthMostRecent(isIn,2);
        
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        
        doSize = zeros(1,3);
        for room = 1:1
            m1 = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]));
            doSize = nanmax([doSize; size(m1)]);
        end
        
        comMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 32]);
        totalMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 4]);
        mfr = nan([doSize(3) 8]);
        
        half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
        
        for first_door = 1:4
            for second_door = 1:4
                comMaps(:,:,:,(first_door-1).*4+second_door) = mkTraceMaps(s.processed.p(:,isInRoom),...
                    gT(:,isInRoom),isSecondMostRecent(second_door,isInRoom) & isMostRecent(first_door,isInRoom) & half(1,isInRoom), ...
                    ones(1,2).*nanmax(doSize(1:2)));
                comMaps(:,:,:,(first_door-1).*4+second_door+16) = mkTraceMaps(s.processed.p(:,isInRoom),...
                    gT(:,isInRoom),isSecondMostRecent(second_door,isInRoom) &isMostRecent(first_door,isInRoom) & ~half(1,isInRoom), ...
                    ones(1,2).*nanmax(doSize(1:2)));
            end
        end
        
        rateComp = nan(32,32,length(totalMaps(1,1,:,1)));
        allComp = nan(32,32);
        rotCorrs = nan(32,32,4);
        for ri = 1:32
            for rj = ri:32
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
        crossAllComp = cat(3,crossAllComp,allComp);
    end
    crossAllComp = crossAllComp(1:16,17:32,:);
 
    crossAllComp = cat(3,nanmean(crossAllComp(:,:,1:6),3), ... 
        nanmean(crossAllComp(:,:,7:12),3), ...
        nanmean(crossAllComp(:,:,13:18),3), ...
        nanmean(crossAllComp(:,:,19:24),3), ...
        nanmean(crossAllComp(:,:,25:30),3) ...
        );
    
    close all
    figure(1)
    set(gcf,'position',[50 50 425 300])
    imagesc(nanmean(crossAllComp(:,:,1:end),3))
    colormap(circshift([linspace(0,1,256)' ...
        [linspace(0,1,128) ones(1,128)]' ...
        linspace(0,1,256)'],[0 -2]))
    caxis([0.0 0.3])
    colorbar
    axis equal
    axis off
    
    
    
    allSim = crossAllComp;
    mask = false(16,16);
    mask(1:4,1:4) = true;
    mask(5:8,5:8) = true;
    mask(9:12,9:12) = true;
    mask(13:16,13:16) = true;
    v = nan(length(allSim(1,1,:)),4);
    for i = 0:3
        v(:,i+1) = help_getMaskedVals(allSim,circshift(mask,[0 i.*4]));
    end
    
    mask = false(16,16);
    mask(1:17:end) = true;
    mask(5:17:end) = true;
    mask(9:17:end) = true;
    mask(13:17:end) = true;
    v2 = nan(length(allSim(1,1,:)),4);
    for i = 0:3
        v2(:,i+1) = help_getMaskedVals(allSim,circshift(mask,[0 i]));
    end
    
    
    figure(2)
    set(gcf,'position',[50 450 500 200])
    subplot(1,2,1)
    mkGraph(v)
    subplot(1,2,2)
    mkGraph(v2)
    
%     figure(3)
%     set(gcf,'position',[50 450 900 200])
%     for i = 1:4
%         tmp = crossAllComp((i-1).*4+1:i.*4,(i-1).*4+1:i.*4);
%         mask = false(4);
%         mask(1:5:end) = true;
%         for q = 0:3
%             tv(:,q+1) = help_getMaskedVals(tmp,circshift(mask,[0 q]));
%         end
%         subplot(1,4,i)
%         mkGraph(tv)
%     end
%     
    root = 'Plots/Summary/FourDoor_2Back/';
    figure(1)
    saveFig(gcf,[root 'RDM'],'pdf');
    figure(2)
    saveFig(gcf,[root 'Similarity'],'pdf');
    
    outP = ['Stats_FourDoor_2Back.txt'];
    fid = fopen(outP,'w');
    fprintf(fid,'\t\t\tSIMILARITY 1-Back\n');
    for i = 1:4
        for j = i+1:4
            [h p ci tstat] = ttest(v(:,i),v(:,j));
            fprintf(fid,['\n' num2str(i) ' to ' num2str(j) ':  ']);
            fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
        end
    end
    
    fprintf(fid,'\n\t\t\tSIMILARITY 2-Back\n');
    for i = 1:4
        for j = i+1:4
            [h p ci tstat] = ttest(v2(:,i),v2(:,j));
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









































