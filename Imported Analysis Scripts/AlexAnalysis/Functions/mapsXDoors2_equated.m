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
        
        [isIn isMostRecent] = isInROI(s.processed.p,s.processed.roi.door);
%         
%         com = isMostRecent(1,:)|isMostRecent(3,:)|isMostRecent(5,:);
%         uni = isMostRecent(2,:)|isMostRecent(4,:)|isMostRecent(6,:);
%         
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        
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
        
        allComp = nan(4,4);
        for doori = 1:4
            for doorj = 1:4
                [a sa] = mkTraceMaps(s.processed.p(:,isInRoom),...
                    gT(:,isInRoom),isMostRecent(doori,isInRoom) & half(1,isInRoom), ...
                    ones(1,2).*nanmax(doSize(1:2)));
                [b sb] = mkTraceMaps(s.processed.p(:,isInRoom),...
                    gT(:,isInRoom),isMostRecent(doorj,isInRoom) & ~half(1,isInRoom), ...
                    ones(1,2).*nanmax(doSize(1:2)));
                
                samp = nanmin(sa,sb);
                [a sa] = mkTraceMaps(s.processed.p(:,isInRoom),...
                    gT(:,isInRoom),isMostRecent(doori,isInRoom) & half(1,isInRoom), ...
                    ones(1,2).*nanmax(doSize(1:2)),samp);
                [b sb] = mkTraceMaps(s.processed.p(:,isInRoom),...
                    gT(:,isInRoom),isMostRecent(doorj,isInRoom) & ~half(1,isInRoom), ...
                    ones(1,2).*nanmax(doSize(1:2)),samp);
                if ~any(~isnan(a)&~isnan(b))
                    continue
                end
                allComp(doori,doorj) = corr(a(~isnan(a)&~isnan(b)),b(~isnan(a)&~isnan(b)));
            end
        end
        
        crossAllComp = cat(3,crossAllComp,allComp);
    end

    crossAllComp = cat(3,nanmean(crossAllComp(:,:,1:6),3), ... 
        nanmean(crossAllComp(:,:,7:12),3), ...
        nanmean(crossAllComp(:,:,13:18),3), ...
        nanmean(crossAllComp(:,:,19:24),3));
    
    close all
    figure(1)
    set(gcf,'position',[50 50 425 300])
    imagesc(nanmean(crossAllComp(:,:,1:end),3))
    colormap(circshift([linspace(0,1,256)' ...
        [linspace(0,1,128) ones(1,128)]' ...
        linspace(0,1,256)'],[0 -2]))
    caxis([0.10 0.4])
    colorbar
    axis equal
    axis off
    
    
    
    allSim = crossAllComp;
    mask = false(4,4);
    mask(1:5:end) = true;
    v = nan(length(allSim(1,1,:)),4);
    for i = 0:3
        v(:,i+1) = help_getMaskedVals(allSim,circshift(mask,[0 i]));
    end
    
    figure(2)
    set(gcf,'position',[50 50 250 200])
    mkGraph(v)
    
    
    
    root = 'Plots/Summary/FourDoor/';
    figure(1)
    saveFig(gcf,[root 'RDM'],'pdf');
    figure(2)
    saveFig(gcf,[root 'Similarity'],'pdf');
    
    outP = ['Stats_FourDoor.txt'];
    fid = fopen(outP,'w');
    fprintf(fid,'\t\t\tSIMILARITY\n');
    for i = 1:4
        for j = i+1:4
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









































