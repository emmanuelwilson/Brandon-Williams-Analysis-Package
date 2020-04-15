function mapsXDoors3(paths)
    crossAllComp = [];
    crossRateComp = [];
    allRotCorrs = [];
    allPredictedRDMs = [];
    for p = paths'
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
        gT = s.processed.trace(s.processed.splithalf.p<=0.05,:);
%         gT = s.processed.trace;
        
        fprintf(['\n\t' p{1} ':  ' num2str(length(gT(:,1)))])

        [isIn isMostRecent indexSinceIn distanceSinceIn] = ...
            isInROI(s.processed.p,s.processed.roi.door(:,[1 2 3]));
        
%         [blah isExit] = isInROI(fliplr(s.processed.p),s.processed.roi.door(:,[1 2 3]));
%         isExit = fliplr(isExit);
        
        [isInRoom blah indexSinceInRoom distanceSinceInRoom] = isInROI(s.processed.p,s.processed.roi.room);

        half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
        
        allComp = nan(6,6);
        
        
%         allMasks = repmat({[]},[1 6]);
%         for i = 1:3
%             allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom)];
%             allMasks{i+3} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom)];
%         end
%         
%         [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks);
        
        for doorA = 1:3
            for doorB = 1:3
                for halfA = 0:1
                    for halfB = 0:1
                        if (halfA).*3+doorA >(halfB).*3+doorB || halfA == halfB
                            continue
                        end
                        [map samp val] = getMatchedSamplingMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),...
                            [isMostRecent(doorA,isInRoom) & halfA==half(1,isInRoom)],...
                            [isMostRecent(doorB,isInRoom) & halfB==half(1,isInRoom)]);
                        allComp((halfA).*3+doorA,(halfB).*3+doorB) = val;
                    end
                end
            end
        end

%         doK = [8 8];
%       
%         doSize = zeros(1,3);
%         for room = 1:1
%             m1 = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
%                 gT(:,[isInRoom(room,:)]));
%             doSize = nanmax([doSize; size(m1)]);
%         end
%         
%         totalMaps = nan([ones(1,2).*nanmax(doSize(1:2)) doSize(3) 4]);
%         for door = 1:3
%             totalMaps(:,:,:,door) = mkTraceMaps(s.processed.p(:,isInRoom),...
%                 gT(:,isInRoom),isMostRecent(door,isInRoom), ...
%                 ones(1,2).*nanmax(doSize(1:2)));
%         end
%         for part = 0:floor(length(totalMaps(1,1,:,1))/prod(doK))
%             
%             figure(1)
%             set(gcf,'position',[50 50 1350 900])
%             for k = 1:prod(doK)
%                 if part.*prod(doK)+k > length(gT(:,1))
%                     break
%                 end
%                 tmp = [totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,2) ...
%                     nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,3)];
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
        crossAllComp = cat(3,crossAllComp,allComp(1:3,4:6));
    end
    
    
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
    caxis([0.0 0.45])
    colorbar
    axis equal
    axis off
    
    
    
    allSim = crossAllComp;
    
    v = nan(length(allSim(1,1,:)),4);
    mask = false(3,3);
    mask(1:4:end) = true;
%     mask(1) = true;
    v(:,1) = help_getMaskedVals(allSim,mask);
    
%     mask = false(3,3);
%     mask(5) = true;
%     v(:,2) = help_getMaskedVals(allSim,mask);
%     
%     
%     mask = false(3,3);
%     mask(9) = true;
%     v(:,3) = help_getMaskedVals(allSim,mask);
    
    mask = false(3,3);
    mask(1,2) = true;
    mask(2,1) = true;
    v(:,2) = help_getMaskedVals(allSim,mask);
    
    mask = false(3,3);
    mask(1,3) = true;
    mask(3,1) = true;
    v(:,3) = help_getMaskedVals(allSim,mask);
    
    mask = false(3,3);
    mask(2,3) = true;
    mask(3,2) = true;
    v(:,4) = help_getMaskedVals(allSim,mask);
    
    figure(2)
    set(gcf,'position',[50 450 250 200])
    mkGraph(v)
    
    
    
    root = 'Plots/Summary/ShortLong/';
    figure(1)
    saveFig(gcf,[root 'RDM'],'pdf');
    figure(2)
    saveFig(gcf,[root 'Similarity'],'pdf');
    
    outP = ['Stats_ShortLong.txt'];
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









































