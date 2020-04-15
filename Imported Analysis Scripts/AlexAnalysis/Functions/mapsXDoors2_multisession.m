function mapsXDoors2(paths)

    
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
    
    crossAllComp = [];
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}) '\n'])
        
        isM = find(ismember(piece,upiece(mi)));
        s = load(paths{isM(1)});
        alignMap = s.alignment.alignmentMap;
        pval = s.alignment.pval;
        mComps = [];
        for di = 1:1:length(isM);
            for dj = di+1:1:length(isM);
                
                fprintf(['\n\t\t( ' num2str(di) ', ' num2str(dj) '; rotation:  ' ...
                    num2str(s.alignment.rotation(di,dj)) ')'])
                
                thresh = 0.05;
                if isempty(alignMap{di,dj}) || nansum(pval{di,dj}<=thresh)<1
                    continue
                end
                tic
                s1 = load(paths{isM(di)});
                s2 = load(paths{isM(dj)});
                
                isInRoom1 = isInROI(s1.processed.p,s1.processed.roi.room);
                isInRoom2 = isInROI(s2.processed.p,s2.processed.roi.room);
                
                [isIn mrd1] = isInROI(s1.processed.p,(s1.processed.roi.door));
                [isIn mrd2] = isInROI(s2.processed.p,(s2.processed.roi.door));
                
                fprintf(['\tCells:  ' num2str(nansum(pval{di,dj}<=thresh)) ...
                    ' (' num2str(100.*nanmean(pval{di,dj}<=thresh)) '%%)']);
                
                P1 = s1.processed.p(:,isInRoom1);
                P2 = s2.processed.p(:,isInRoom2);
                T1 = s1.processed.trace(alignMap{di,dj}(pval{di,dj}<=thresh,1),isInRoom1);
                T2 = s2.processed.trace(alignMap{di,dj}(pval{di,dj}<=thresh,2),isInRoom2);
                T = [T1 T2];
                P = [P1 P2];
                
                allMasks = repmat({[]},[1 8]);
                for i = 1:4
                    allMasks{i} = [mrd1(i,isInRoom1) false(1,length(mrd2(i,isInRoom2)))];
                    allMasks{i+4} = [false(1,length(mrd1(i,isInRoom1))) mrd2(i,isInRoom2)];
                end

                [map samp allComp ivals] = getMatchedMapsNMasks(P,T,allMasks);
                mComps = cat(3,mComps,ivals(1:4,5:8,:));
                
                totalMaps = [];
                for i = 1:4
                    if isempty(totalMaps)
                        totalMaps = getMatchedMapsNMasks(P,T,{[mrd1(i,isInRoom1) mrd2(i,isInRoom2)]});
                    else
                        totalMaps = cat(2,totalMaps,nan([length(totalMaps(:,1,1)) 1 length(totalMaps(1,1,:))]),...
                            getMatchedMapsNMasks(P,T,{[mrd1(i,isInRoom1) mrd2(i,isInRoom2)]}));
                    end
                end
                
                doK = [8 4];
        
                for part = 0:floor(length(totalMaps(1,1,:,1))/prod(doK))

                    figure(1)
                    set(gcf,'position',[50 50 1350 900])
                    for k = 1:prod(doK)
                        if part.*prod(doK)+k > length(T(:,1))
                            break
                        end

                        subplot(doK(1),doK(2),k)    
                        imagesc(totalMaps(:,:,part.*prod(doK)+k))
                        colormap jet
                        alpha(double(~isnan(totalMaps(:,:,part.*prod(doK)+k))))
                        axis equal
                        axis off    
                    end

                    slashInds = find(ismember(upiece{mi},'/'));
                    outP = ['Plots/DifferentiatedMaps/' upiece{mi}(slashInds(1)+1:end) '/' spiece{isM(di)} ...
                        '_' spiece{isM(dj)} '_Partition_' num2str(part+1)];
                    saveFig(gcf,outP,'tiff')
                    saveFig(gcf,outP,'pdf')
                    close all
                    drawnow
                end
            end
        end
        crossAllComp = cat(3,crossAllComp,nanmean(mComps,3));
    end

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
    
    
    
    allSim = mComps;
    mask = false(4,4);
    mask(1:5:end) = true;
    v = nan(length(allSim(1,1,:)),4);
    for i = 0:3
        v(:,i+1) = help_getMaskedVals(allSim,circshift(mask,[0 i]));
    end
    
    figure(2)
    set(gcf,'position',[50 450 500 200])
%     v2 = load('predict','v');
    mkGraph([{v(:,1)} {v(:,2)} {v(:,3)} {v(:,4)}])
    set(gca,'ylim',[0 0.4])
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









































