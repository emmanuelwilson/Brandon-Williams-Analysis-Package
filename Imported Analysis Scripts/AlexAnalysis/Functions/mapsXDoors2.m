function mapsXDoors2(paths)
    clc
    crossAllComp = [];
    crossRateComp = [];
    allRotCorrs = [];
    for p = paths'
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
%         gT = s.processed.trace(s.processed.splithalf.p<=0.01,:);
         gT = s.processed.trace(s.processed.splithalf.within.p<=0.05,:);
%         gT = s.processed.trace;

        fprintf(['\n\t' p{1} ':  ' num2str(length(gT(:,1)))])

        [isIn isMostRecent] = isInROI(s.processed.p,(s.processed.roi.door));
%         isMostRecent = fliplr(isMostRecent);
%         isMostRecent = nthMostRecent(isIn,2);
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        
        half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
        
        allMasks = repmat({[]},[1 8]);
        for i = 1:4
            allMasks{i} = [isMostRecent(i,isInRoom) & half(:,isInRoom)];
            allMasks{i+4} = [isMostRecent(i,isInRoom) & ~half(:,isInRoom)];
        end
        
        [maps samp allComp ivals] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks);
        
%         allComp = nan(8,8);
%         for doorA = 1:4
%             for doorB = 1:4
%                 for halfA = 0:1
%                     for halfB = 0:1
%                         if (halfA).*4+doorA > (halfB).*4+doorB || halfA == halfB
%                             continue
%                         end
%                         [map samp val] = getMatchedSamplingMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),...
%                             [isMostRecent(doorA,isInRoom) & halfA==half(1,isInRoom)],...
%                             [isMostRecent(doorB,isInRoom) & halfB==half(1,isInRoom)]);
%                         allComp((halfA).*4+doorA,(halfB).*4+doorB) = val;
%                     end
%                 end
%             end
%         end
        
%         doSize = zeros(1,3);
%         for room = 1:1
%             m1 = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
%                 gT(:,[isInRoom(room,:)]));
%             doSize = nanmax([doSize; size(m1)]);
%         end
%         
%         comMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 8]);
%         totalMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 4]);
%         mfr = nan([doSize(3) 8]);
%         
%         half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
%         
%         for door = 1:4
%             comMaps(:,:,:,door) = mkTraceMaps(s.processed.p(:,isInRoom),...
%                 gT(:,isInRoom),isMostRecent(door,isInRoom) & half(1,isInRoom), ...
%                 ones(1,2).*nanmax(doSize(1:2)));
%             comMaps(:,:,:,door+4) = mkTraceMaps(s.processed.p(:,isInRoom),...
%                 gT(:,isInRoom),isMostRecent(door,isInRoom) & ~half(1,isInRoom), ...
%                 ones(1,2).*nanmax(doSize(1:2)));
%             
%             mfr(:,door) = nanmean(gT(:,isMostRecent(door,:) & isInRoom & half),2);
%             mfr(:,door+4) = nanmean(gT(:,isMostRecent(door,:) & isInRoom & ~half),2);
% 
%             
%             totalMaps(:,:,:,door) = mkTraceMaps(s.processed.p(:,isInRoom),...
%                 gT(:,isInRoom),isMostRecent(door,isInRoom), ...
%                 ones(1,2).*nanmax(doSize(1:2)));
%         end
%        
%         
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
%                 tmp = [[totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,2)] ; ...
%                     nan(1,size(totalMaps,2).*2 + 1); ...
%                     [totalMaps(:,:,part.*prod(doK)+k,3) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,4)]];
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
        
%         crossRateComp = cat(3,crossRateComp,rateComp);
        crossAllComp = cat(3,crossAllComp,ivals(1:4,5:8,:));
%         allRotCorrs = cat(4,allRotCorrs,bmr);
    end
    
    
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
    %% Average within animal
% % %     piece = [];
% % %     for i = 1:length(paths)
% % %         ind = find(ismember(paths{i},'/'),1,'last')-1;
% % %         piece = [piece; {paths{i}(1:ind)}];
% % %     end
% % %     upiece = unique(piece);
% % %     val = [];
% % %     for i = 1:length(upiece)
% % %         val = cat(3,val,nanmean(crossAllComp(:,:,ismember(piece,upiece(i))),3));
% % %     end
% % %     crossAllComp = val;
    
    close all
    figure(1)
    set(gcf,'position',[50 50 425 300])
    imagesc(nanmean(crossAllComp(:,:,1:end),3))
    colormap(circshift([linspace(0,1,256)' ...
        [linspace(0,1,128) ones(1,128)]' ...
        linspace(0,1,256)'],[0 -2]))
    caxis([0.0 0.45])
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
    set(gcf,'position',[50 450 500 200])
    v2 = load('predict','v');
    mkGraph([{v(:,1)} {v(:,2)} {v(:,3)} {v(:,4)}; ...
        {v2.v(:,1)} {v2.v(:,2)} {v2.v(:,3)} {v2.v(:,4)}]')
    set(gca,'ylim',[0 0.45])
%     save('predict','v');
    
%     figure(3)
%     set(gcf,'position',[50 50 800 800])
%     for i = 1:4
%         for j = 1:4
%             if i == j
%                 continue
%             end
%             subplot(4,4,(i-1).*4+j)
%             tv = permute(allRotCorrs(i,j,:),[3 1 2]);
%             mkGraph([tv==1 tv==2 tv==3 tv==4],[0 90 180 270])
%             set(gca,'ylim',[0 1])
%         end
%     end

    root = 'Plots/Summary/FourDoor/';
    figure(1)
    saveFig(gcf,[root 'RDM'],'pdf');
    figure(2)
    saveFig(gcf,[root 'Similarity'],'pdf');
%     figure(3)
%     saveFig(gcf,[root 'BestMatchRotations'],'pdf');
%     saveFig(gcf,[root 'BestMatchRotations'],'tiff');
    
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









































