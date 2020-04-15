function mapsXDoors5_parts(paths)
    crossAllComp = [];
    crossDataAmount = [];
    labels = [];
    groups = [{'Saline'} {'CNO'}];
    for pi = 1:length(paths)
        p = paths(pi);
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
%         gT = s.processed.trace(s.processed.splithalf.p<=0.05,:);
        gT = s.processed.trace(s.processed.splithalf.within.p<=0.05,:);
%         gT = s.processed.trace;
        
%         if length(gT(:,1))<10
%             paths(pi) = [];
%             continue
%         end
        
        labels = [labels; p {length(gT(:,1))}];
        
        fprintf(['\n\t' p{1} ':  ' num2str(length(gT(:,1)))])

        [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        a = nanmin(indexSinceIn).*(1./30);
        oIndexSinceIn = indexSinceIn;
%         doRange = [0:5:30];
        doRange = [0:1:10];
        rangeComp = nan(2,length(doRange));
        dataAmount = nan(1,length(doRange));
        for curR = 1:length(doRange)
            indexSinceIn = oIndexSinceIn;
%             indexSinceIn = nanmin(indexSinceIn).*(1./30) >= doRange(curR) & ...
%                 nanmin(indexSinceIn).*(1./30) < doRange(curR+1);
            indexSinceIn = nanmin(indexSinceIn).*(1./30) >= doRange(curR);
            dataAmount(curR) = nansum(indexSinceIn);
            if ~any(indexSinceIn)
                continue
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

            allMasks = repmat({[]},[1 6]);
            for i = 1:2
                allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom) & indexSinceIn(1,isInRoom)];
                allMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom) & indexSinceIn(1,isInRoom)];
            end

            [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks);

            
%             allComp = nan(4,4);
%             for doorA = 1:2
%                 for doorB = 1:2
%                     for halfA = 0:1
%                         for halfB = 0:1
%                             if (halfA).*2+doorA > (halfB).*2+doorB || halfA == halfB
%                                 continue
%                             end
%                             [map samp val] = getMatchedSamplingMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),...
%                                 [isMostRecent(doorA,isInRoom) & halfA==half(1,isInRoom) & indexSinceIn(1,isInRoom)],...
%                                 [isMostRecent(doorB,isInRoom) & halfB==half(1,isInRoom) & indexSinceIn(1,isInRoom)]);
%                             allComp((halfA).*2+doorA,(halfB).*2+doorB) = val;
%                         end
%                     end
%                 end
%             end

            v = nan(1,2);
            mask = false(2,2);
            mask(1:3:end) = true;
            v(1,1) = help_getMaskedVals(allComp(1:2,3:4),mask);

            mask = false(2,2);
            mask(1,2) = true;
            mask(2,1) = true;
            v(1,2) = help_getMaskedVals(allComp(1:2,3:4),mask);
    
            rangeComp(:,curR) = v';
        end
        crossAllComp = cat(3,crossAllComp,rangeComp);
        crossDataAmount = [crossDataAmount; dataAmount];
    end
%     crossAllComp = crossAllComp(:,:,end:-1:1);

%     close all
%     figure(1)
%     set(gcf,'position',[50 50 425 300])
%     imagesc(nanmean(crossAllComp(:,:,1:end),3))
%     colormap(circshift([linspace(0,1,256)' ...
%         [linspace(0,1,128) ones(1,128)]' ...
%         linspace(0,1,256)'],[0 -2]))
%     caxis([0.0 0.5])
%     colorbar
%     axis equal
%     axis off

    gm =  nan(length(paths),1);
    for i = 1:length(paths)
        gm(i) = find(ismember(groups,paths{i}(find(ismember(paths{i},'_'),1,'last')+1:end-4)));
    end

    if all(gm==0)

        piece = [];
        for i = 1:length(paths)
            ind = find(ismember(paths{i},'/'),1,'last')-1;
            piece = [piece; {paths{i}(1:ind)}];
        end
        upiece = unique(piece);
        val = [];
        byAnimal = [];
        bySession = nan(length(upiece),21);
        for i = 1:length(upiece)
            val = cat(3,val,nanmean(crossAllComp(:,:,ismember(piece,upiece(i))),3));
            byAnimal{i} = crossAllComp(:,:,ismember(piece,upiece(i)));
            bySession(i,1:nansum(ismember(piece,upiece(i)))) = ...
                permute(diff(byAnimal{i}(:,1,:),[],1),[1 3 2]);
        end
        step = 3;
        doPlot = [];
        for i = 1:step:length(bySession(1,:))
            doPlot = [doPlot nanmean(bySession(:,i:min(i+step-1,length(bySession(1,:)))),2)];
        end
        figure(4)
        set(gcf,'position',[50 50 600 600])
        subplot(2,1,1)
        mkGraph(doPlot)
        subplot(2,1,2)
        mkLine({bySession})

    %     crossAllComp = val;
    %     [labels crossAllComp]

        figure(2)
        set(gcf,'position',[50 450 450 450])
        subplot(3,1,[1:2])
        tmp = permute(crossAllComp,[3 2 1]);
        tmp2 = permute(diff(crossAllComp,[],1),[3 2 1]);
        tmp = tmp(:,1:end,:);
        tmp2 = -tmp2(:,1:end);
        h = mkLine([{tmp(:,:,1)} {tmp(:,:,2)} {tmp2}],doRange);
        legend(h,[{'Same entry'} {'Other entry'} {'Difference'}],'location','northeast')
        set(gca,'ylim',[-0.2 0.55])
        tmp3 = bsxfun(@rdivide,crossDataAmount,crossDataAmount(:,1));
        subplot(3,1,3)
        mkLine({tmp3(:,1:end)},doRange(1:end));
        hold on
        plot([doRange(1) doRange(end)],[1 1],'color','k','linestyle','-',...
            'linewidth',1)
        set(gca,'ylim',[0 1.1])
        drawnow

        root = 'Plots/Summary/TwoBig/';
        figure(2)
        saveFig(gcf,[root 'RemappingWithinTrial'],'pdf');
        saveFig(gcf,[root 'RemappingWithinTrial'],'tiff');

        figure(4)
        saveFig(gcf,[root 'RemappingAcrossTrial'],'pdf');
        saveFig(gcf,[root 'RemappingAcrossTrial'],'tiff');

        [h p ci tstat] = ttest(tmp2,0)

        %%%%%%%%%%%%%%%%%%% By groups for CNO
        
    else
    
        tmp = permute(crossAllComp,[3 2 1]);
        tmp2 = permute(diff(crossAllComp,[],1),[3 2 1]);
        tmp = tmp(:,1:end,:);
        tmp2 = -tmp2(:,1:end);

        figure(5)
        set(gcf,'position',[50 50 350 250])
        h = mkLine([{tmp2(gm==1,:)} {tmp2(gm==2,:)}],doRange);
        set(gca,'xlim',[-0.5 doRange(end)+0.5])


        figure(6)
        for group = 1:length(groups)
            set(gcf,'position',[50 450 700 450])
            subplot(3,2,[1 3]+[group-1])
            tmp = permute(crossAllComp,[3 2 1]);
            tmp2 = permute(diff(crossAllComp,[],1),[3 2 1]);
            tmp = tmp(:,1:end,:);
            tmp2 = tmp2(:,1:end);
            h = mkLine([{tmp(gm==group,:,1)} {tmp(gm==group,:,2)}],doRange);
    %         legend(h,[{'Same entry'} {'Other entry'} {'Difference'}],'location','northeast')
            set(gca,'ylim',[-0.0 0.65])
            set(gca,'xlim',[-0.5 doRange(end)+0.5])
            tmp3 = bsxfun(@rdivide,crossDataAmount,crossDataAmount(:,1));
            subplot(3,2,5+group-1)
            mkLine({tmp3(:,1:end)},doRange(1:end));
            hold on
            plot([doRange(1) doRange(end)],[1 1],'color','k','linestyle','-',...
                'linewidth',1)
            set(gca,'ylim',[0 1.1])
            set(gca,'xlim',[-0.5 doRange(end)+0.5])
            drawnow
        end


        root = 'Plots/Summary/DREADDs_TwoSmall/';
        figure(5)
        saveFig(gcf,[root 'CNO_RemappingAcrossTime_Diff'],'pdf');
        saveFig(gcf,[root 'CNO_RemappingAcrossTime_Diff'],'tiff');
        figure(6)
        saveFig(gcf,[root 'CNO_RemappingAcrossTime'],'pdf');
        saveFig(gcf,[root 'CNO_RemappingAcrossTime'],'tiff');

        outP = ['Stats_DREADDs_TwoSmall.txt'];
        fid = fopen(outP,'w');
        fprintf(fid,'\t\t\tSIMILARITY ACROSS TIME\n');

        [h p ci tstat] = ttest(tmp2(gm==1,:),tmp2(gm==2,:));
        fprintf(fid,['\n\nSaline vs. CNO']);
        fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',[tstat.df; tstat.tstat; p]);
        [h p ci tstat] = ttest(tmp2(gm==1,:),0);
        fprintf(fid,['\n\nSaline vs. 0']);
        fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',[tstat.df; tstat.tstat; p]);
        [h p ci tstat] = ttest(tmp2(gm==2,:),0);
        fprintf(fid,['\n\nCNO vs. 0']);
        fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',[tstat.df; tstat.tstat; p]);

        fclose all
    end
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end









































