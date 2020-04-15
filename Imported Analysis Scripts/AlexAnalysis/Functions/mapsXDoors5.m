function mapsXDoors5(paths)
    warning off all
    
    
    labels = [{'Saline'} {'CNO'} {'No Injection'}];
    afw = repmat({[]},[length(labels) 4]);
    apfr = repmat({[]},[length(labels) 4]);
    amfr = repmat({[]},[1 length(labels)]);
    icorr = repmat({[]},[length(labels) 4]);
    cellCounts = repmat({[]},[1 length(labels)]);
    crossAllComp = repmat({[]},[1 length(labels)]);
    allIValDiffs = repmat({[]},[length(paths) 1]);
    allSNR = [];
    for p = paths'
        s = load(p{1});
        group = find(ismember(labels,p{1}(find(ismember(p{1},'_'),1,'last')+1:end-4)));
        if isempty(group)
            group = 3;
        end
        vel = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
%         gT = s.processed.trace(s.processed.splithalf.room.p<=0.01,:);
% % %         snr = nanmean(s.processed.trace,2)./nanstd(s.processed.trace,[],2);
        
%         allSNR = [allSNR; snr s.processed.splithalf.wholemap.p];

% % %         [a b c] = mkTraceMaps(s.processed.p,s.processed.trace);
% % %         b = b./nansum(b(:));
% % %         si = permute(nansum(nansum(bsxfun(@times,b,(c./repmat(nanmean(nanmean(c,1),2),[size(b)])) .* ...
% % %             log(c./repmat(nanmean(nanmean(c,1),2),[size(b)]))),1),2),[3 1 2]);
% % %         si = s.processed.splithalf.wholemap_si.p;
        % s.processed.splithalf.roomXdoors.p<=0.05 & si > 
%         inds = (s.processed.splithalf.wholemap_si.p<=0.05 & s.processed.splithalf.wholemap.p<=0.05);% & ~s.processed.exclude.SFPs);
        cellType = nan(1,4);
        if isfield(s.processed,'exclude')
            inds = (s.processed.splithalf.roomXdoors.p<=0.05 & s.processed.exclude.SFPs);
            inc = s.processed.exclude.SFPs;
            cellType = [nansum(inc & ~s.processed.splithalf.roomXdoors.p<=0.05 & ~s.processed.splithalf.hallwayXdoors.p<=0.05) ...
                nansum(inc & s.processed.splithalf.roomXdoors.p<=0.05 & ~s.processed.splithalf.hallwayXdoors.p<=0.05) ...
                nansum(inc & ~s.processed.splithalf.roomXdoors.p<=0.05 & s.processed.splithalf.hallwayXdoors.p<=0.05) ...
                nansum(inc & s.processed.splithalf.roomXdoors.p<=0.05 & s.processed.splithalf.hallwayXdoors.p<=0.05)];
        else
            inds = (s.processed.splithalf.roomXdoors.p<=0.05);% & ~s.processed.exclude.SFPs);
            
            cellType = [nansum(s.processed.splithalf.roomXdoors.p>0.05 & s.processed.splithalf.hallwayXdoors.p>0.05) ...
                nansum(s.processed.splithalf.roomXdoors.p<=0.05 & s.processed.splithalf.hallwayXdoors.p>0.05) ...
                nansum(s.processed.splithalf.roomXdoors.p>0.05 & s.processed.splithalf.hallwayXdoors.p<=0.05) ...
                nansum(s.processed.splithalf.roomXdoors.p<=0.05 & s.processed.splithalf.hallwayXdoors.p<=0.05)];
        end
        
        cellCounts{group} = [cellCounts{group}; cellType];
        
        gT = s.processed.trace(inds,:);
% %         [a b] = sort(1-s.processed.splithalf.rate.room.p(inds,1),'descend');
% %         gT = gT(b,:);
% % %         snr = nanmean(gT,2)./nanstd(gT,[],2);
%         gT(snr>0.20,:) = [];
%         gT = s.processed.trace(:,:);
%         gT = gT>0;


% %              num2str(100.*nanmean(s.processed.splithalf.roomXdoors.p(s.processed.exclude.SFPs)<=0.05))

         fprintf(['\n\t' p{1} ':  ' num2str(length(gT(:,1)))])

        [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        a = nanmin(indexSinceIn).*(1./30);
%         thresh = median(a(isInRoom));
        thresh = -1;
%         thresh = 5;
        
        indexSinceIn = nanmin(indexSinceIn).*(1./30) >= thresh;
        
        half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;

        allMasks = repmat({[]},[1 4]);
        for i = 1:2
            allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom) & vel(1,isInRoom)>-2];
            allMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom) & vel(1,isInRoom)>-2];
%             h2 = cumsum(isMostRecent(i,isInRoom));
%             allMasks{i} = [isMostRecent(i,isInRoom) & h2<h2(end)./2];
%             allMasks{i+2} = [isMostRecent(i,isInRoom) & ~(h2<h2(end)./2)];
        end
        
        [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks);
%         [allComp ival mfr] = getMatchedSamplingValues(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks);
        
        v = help_getMaskedVals(allComp(1:2,3:4,:),[true false; false true]) - ...
            help_getMaskedVals(allComp(1:2,3:4,:),[false true; true false]);
        
        s.processed.similarity.roomXdoor.pop = v;
        
%         v = help_getMaskedVals(mfr(1:2,3:4,:),[true false; false true]) - ...
%             help_getMaskedVals(mfr(1:2,3:4,:),[false true; true false]);
        s.processed.similarity.roomXdoor.ival = v;
        
        v = help_getMaskedVals(allComp(1:2,3:4,:),[true false; false true]) - ...
            help_getMaskedVals(allComp(1:2,3:4,:),[false true; true false]);
        
        [blah1 ivals mfr pfr] = getMatchedSamplingValues(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks);
        apfr(group,:) = [{[apfr{group,1}; help_getMaskedVals(mfr(1:2,3:4,:),[true false; false true])]} ...
            {[]} {[apfr{group,3}; help_getMaskedVals(mfr(1:2,3:4,:),[false true; true false])]} {[]}];
        icorr(group,:) = [{[icorr{group,1}; help_getMaskedVals(ivals(1:2,3:4,:),[true false; false true])]} ...
            {[]} {[icorr{group,3}; help_getMaskedVals(ivals(1:2,3:4,:),[false true; true false])]} {[]}];
        
%         crossRateComp = cat(3,crossRateComp,rateComp);
        crossAllComp{group} = cat(3,crossAllComp{group},allComp(1:2,3:4));
%         crossAllComp{group} = cat(3,crossAllComp{group},nanmean(ivals(1:2,3:4,:),3));
        

        [w a] = help_getfieldshifts(map);
        nmap = cat(4,map(:,:,randperm(size(map,3)),1),map(:,:,randperm(size(map,3)),2),...
            map(:,:,randperm(size(map,3)),3),map(:,:,randperm(size(map,3)),4));
        [nw na] = help_getfieldshifts(nmap);
        afw(group,:) = [{[afw{group,1}; w]} {[afw{group,2}; nw]} ...
            {[afw{group,3}; a]} [afw{group,4}; na]];
      
%         h1 = cumsum(isMostRecent(1,isInRoom))<nansum(isMostRecent(1,isInRoom))./2;
%         h2 = cumsum(isMostRecent(2,isInRoom))<nansum(isMostRecent(2,isInRoom))./2;
%         
% % %         tgt = gT(:,isInRoom);
% % %         tmp = [nanmean(tgt(:,isMostRecent(1,isInRoom)&h1),2) ... 
% % %             nanmean(tgt(:,isMostRecent(2,isInRoom)&h2),2) ...
% % %             nanmean(tgt(:,isMostRecent(1,isInRoom)&~h1),2) ...
% % %             nanmean(tgt(:,isMostRecent(2,isInRoom)&~h2),2)];

% % %         tmp = [nanmean(gT(:,isMostRecent(1,:)&isInRoom&half),2) ... 
% % %             nanmean(gT(:,isMostRecent(2,:)&isInRoom&half),2) ...
% % %             nanmean(gT(:,isMostRecent(1,:)&isInRoom&~half),2) ...
% % %             nanmean(gT(:,isMostRecent(2,:)&isInRoom&~half),2)];
% % %         
% % %         tmp(:,1:2) = bsxfun(@rdivide,tmp(:,1:2),nanmax(tmp(:,1:2),[],2));
% % %         tmp(:,3:4) = bsxfun(@rdivide,tmp(:,3:4),nanmax(tmp(:,3:4),[],2));
% % % 
% % %         tmpC = [nanmean(abs([diff(tmp(:,[1 3]),[],2) diff(tmp(:,[2 4]),[],2)]),2) ...
% % %             nanmean(abs([diff(tmp(:,[1 4]),[],2) diff(tmp(:,[2 3]),[],2)]),2)];
% % %         
% % %         apfr(group,[1 3]) = [{[apfr{group,1}; nanmean(abs([diff(tmp(:,[1 3]),[],2) diff(tmp(:,[2 4]),[],2)]),2)]} ...
% % %             {[apfr{group,3}; nanmean(abs([diff(tmp(:,[2 3]),[],2) diff(tmp(:,[1 4]),[],2)]),2)]}];
        
        amfr{group} = [amfr{group}; nanmean(gT,2)];
        
%         apfr(group,:) = [{[apfr{group,1}; nanmean([abs(tmp(:,1)-tmp(:,3))./nanmax(tmp(:,1),tmp(:,3)) ...
%             abs(tmp(:,2)-tmp(:,4))./nanmax(tmp(:,2),tmp(:,4))],2)]} ... 
%             {[apfr{group,2}; nanmean([abs(tmpB(:,1)-tmpB(:,3))./nanmax(tmpB(:,1),tmpB(:,3)) ...
%             abs(tmpB(:,2)-tmpB(:,4))./nanmax(tmpB(:,2),tmpB(:,4))],2)]} ... 
%             {[apfr{group,3}; nanmean([abs(tmp(:,1)-tmp(:,4))./nanmax(tmp(:,1),tmp(:,4)) ...
%             abs(tmp(:,2)-tmp(:,3))./nanmax(tmp(:,2),tmp(:,3))],2)]} ...
%             {[apfr{group,4}; nanmean([abs(tmpB(:,1)-tmpB(:,4))./nanmax(tmpB(:,1),tmpB(:,4)) ...
%             abs(tmpB(:,2)-tmpB(:,3))./nanmax(tmpB(:,2),tmpB(:,3))],2)]}];
        
%         [a b] = sort(diff(tmpC,[],2),'descend');
%         tmp = [help_getMaskedVals(mfr(1:2,3:4,:),[false true; true false]) - ...
%             help_getMaskedVals(mfr(1:2,3:4,:),[true false; false true])];
%         [a b] = sort(tmp,'descend');
%         gT = gT(b,:);
     
        doSize = zeros(1,2);
        for room = 1:1
            m1 = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]));
            doSize = nanmax([doSize; size(m1(:,:,1))]);
        end
        
        clear totalMaps
        for door = 1:2
            
            totalMaps(:,:,:,door) = mkTraceMaps(s.processed.p(:,isInRoom),...
                gT(:,isInRoom),isMostRecent(door,isInRoom), ...
                (doSize(1:2)));
        end
        
        totalMaps = map;

        doK = [8 4];
        
        for part = 0:floor(length(totalMaps(1,1,:,1))/prod(doK))
            
            figure(1)
            set(gcf,'position',[50 50 900 1350])
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(gT(:,1))
                    break
                end
%                 tmp = [totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,2)];
                
                tmp = [totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,3); ...
                    nan(1,length(totalMaps(1,:,1,1)).*2+1); ...
                    totalMaps(:,:,part.*prod(doK)+k,2) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,4)];

                
                subplot(doK(1),doK(2),k)
                
% % % %                 step = 3;
% % % %                 plot(s.processed.p(2,1:step:end),s.processed.p(1,1:step:end), ...
% % % %                     'color',[0.0 0.0 0.0],'linewidth',1) %[0.8 0.8 0.8]
% % % %                 hold on
% % % %                 
% % % %                 plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))&isMostRecent(1,:)),...
% % % %                     s.processed.p(1,logical(gT(part.*prod(doK)+k,:))&isMostRecent(1,:)),...
% % % %                     'color',[0.5 0.1 1],'linestyle','none','marker','o',...
% % % %                     'markerfacecolor',[0.5 0.1 1],'markersize',2);
% % % % 
% % % %                 plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))&isMostRecent(2,:)),...
% % % %                     s.processed.p(1,logical(gT(part.*prod(doK)+k,:))&isMostRecent(2,:)),...
% % % %                     'color',[0.1 1 1],'linestyle','none','marker','o',...
% % % %                     'markerfacecolor',[0.1 1 1],'markersize',2);

                imagesc(tmp)
                colormap jet
                caxis([0 nanmax(tmp(:))])
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
        
%         save(p{1},'-struct','s','-v7.3');
    end
    
% % %     vals = [];
% % %     for i = 0:0.02:0.5
% % %         vals = [vals nanmean(allSNR(allSNR(:,1)>i & allSNR(:,1)<i+0.02,2)<0.05)];
% % %     end
% % %     plot(0:0.02:0.5,vals)
    
    close all
    
    labels(cellfun(@isempty,crossAllComp)) = [];
    apfr(cellfun(@isempty,crossAllComp),:) = [];
    icorr(cellfun(@isempty,crossAllComp),:) = [];
    afw(cellfun(@isempty,crossAllComp),:) = [];
    cellCounts(cellfun(@isempty,crossAllComp)) = [];
%     amfr(cellfun(@isempty,crossAllComp),:) = [];
    crossAllComp(cellfun(@isempty,crossAllComp)) = [];
    
%     mkPie(cellCounts)
    
    %% Average within animal
    piece = [];
    ag = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        groupPiece = find(ismember(labels,paths{i}(find(ismember(paths{i},'_'),1,'last')+1:end-4)));
        if isempty(groupPiece)
            groupPiece = 1;
        end
        ag = [ag; groupPiece];
    end
    if ~isempty(ag)
        upiece = unique(piece);
        pairNum = nan(length(piece),2);
        count = zeros(1,length(upiece));
        for i = 1:length(piece)
            count(ismember(upiece,piece(i))) = count(ismember(upiece,piece(i)))+1;
            pairNum(i,ag(i)) = ceil(count(ismember(upiece,piece(i)))./2);
        end
    end

    
%     tmp = cellfun(@nanmedian,allIValDiffs);
%     tmp(isnan(pairNum(:,1))) = [];
%     mkGraph([{tmp(1:4)} {tmp(5:10)} {tmp(11:16)}]);
    
    %%%%%%%%%
    
    av = repmat({[]},[1 length(labels)]);
    toPlot = repmat({[]},[2 length(labels)]);
    for group = 1:length(labels)
    
        allSim = crossAllComp{group};
        if isempty(allSim)
            continue
        end
        
        v = nan(length(allSim(1,1,:)),2);
        mask = false(2,2);
        mask(1:3:end) = true;
        v(:,1) = help_getMaskedVals(allSim,mask);

        mask = false(2,2);
        mask(1,2) = true;
        mask(2,1) = true;
        v(:,2) = help_getMaskedVals(allSim,mask);
        av{group} = v;
        toPlot{1,group} = v(:,1);
        toPlot{2,group} = v(:,2);
    end
    
    figure(2)
    set(gcf,'position',[50 450 400 200])
    subplot(1,3,1:2)
    mkGraph(toPlot,labels,[{'Same Entry'} {'Diff. Entry'}])
%     set(gca,'ylim',[0 0.7])
    ylabel('Correlation (r)')
    title('Split-half map similarity')
    
	subplot(1,3,3)  
    av = cellfun(@fliplr,av,'uniformoutput',false);
    tmp = cellfun(@diff,av,repmat({[]},[1 length(labels)]),repmat({2},[1 length(labels)]),'uniformoutput',false);
    for i = 1:length(tmp)
        tmp{i} = [100.*tmp{i}./av{i}(:,2)];
    end
    mkGraph(tmp(~cellfun(@isempty,tmp)),[])
    set(gca,'ylim',[-5 30])
    ylabel('Correlation Difference (%)')

    figure(4)
    set(gcf,'position',[50 350 350 250])
    animalNormed = repmat({[]},[length(upiece) length(labels)]);

%     animalNormedHallway = repmat({[]},[length(upiece) length(labels)]);
    for i = 1:length(upiece)
        for j = 1:length(labels)
            animalNormed{i,j} = tmp{j}(ismember(piece(ag==j),upiece(i)));

%             animalNormedHallway{i,j} = tmp2.tmp{j}(ismember(piece(ag==j),upiece(i)));
        end
%         if length(animalNormed{i,1})==length(animalNormed{i,2})
%             animalNormed{i,3} = animalNormed{i,1}-animalNormed{i,2};
%         end
    end
%     hist(tmp2,[-0.15:0.015:0.15])
    mkGraph(animalNormed(:,any(~cellfun(@isempty,animalNormed),1))');
    
    figure(21)
    set(gcf,'position',[500 50 800 300])
    subplot(1,2,1)
    cumHist({apfr{1,3}-apfr{1,1}},[-1.5:0.02:1.5])
    
%     root = 'Plots/Summary/DREADDs_TwoSmall/';
    root = 'Plots/Summary/TwoSmall_CA1/';
    for i = 1:length(labels)
        figure(10+i)
        h = cumHist(apfr(i,[1 3]),[0:0.025:1]);
        xlabel('Mean Rate Change (%)')
        ylabel('Count')
        legend([h{1}(1) h{2}(1)],[{'Within'} ...
            {'Across'} {'Shuffled'}],'location','northeast','fontname','arial',...
            'fontsize',9,'fontweight','bold','color','none','box','off')
        set(gcf,'position',[50 400 350 300])
        set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
        title(labels{i})
        drawnow
        saveFig(gcf,[root 'RateRemapping_' labels{i}],[{'tiff'} {'pdf'}]);
        
        figure(20+i)
        h = cumHist(afw(i,[1 3]),[0:1:35]);
%         h = cumHist(icorr(i,[1 3]),[-1:0.025:1]);
        xlabel('Field peak shift (cm)')
        ylabel('Count')
        legend([h{1}(1) h{2}(1)],[{'Within'} ...
            {'Across'} {'Shuffled'}],'location','northeast','fontname','arial',...
            'fontsize',9,'fontweight','bold','color','none','box','off')
        set(gcf,'position',[50 50 350 300])
        drawnow
        title(labels{i})
        saveFig(gcf,[root 'FieldShift_' labels{i}],[{'tiff'} {'pdf'}]);
    end
    
    figure(20)
    set(gcf,'position',[500 50 800 300])
    subplot(1,2,1)
    hist(apfr{1,3}-apfr{1,1},[-1:0.02:1])
    subplot(1,2,2)
    hist(apfr{2,3}-apfr{2,1},[-1:0.02:1])
%     cumHist([{apfr{1,3}-apfr{1,1}} {apfr{2,3}-apfr{2,1}}],[-0.8:0.02:0.8])

%     scatter(apfr{1,1},apfr{1,3})
%     hold on
%     plot([0 1],[0 1],'color','k','linestyle','--')
    
    [h p ci tstat] = ttest(tmp{1},tmp{2});
    fprintf('\n\tPaired t-test:\n\t\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    [h p ci tstat] = ttest2(tmp{1},tmp{2});
    fprintf('\n\t2-Sample t-test:\n\t\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);

    %%%%%%%%%%%%%%%%%%% Versus Hallway
    
% % %     tmp2 = load('DREADDs_TwoSmall_HallwayRemapping');
% % %     
% % %     figure(6)
% % %     set(gcf,'position',[50 750 500 250])
% % %     subplot(1,2,1)
% % %     plot(tmp2.tmp{1},tmp{1},'linestyle','none','marker','o','markersize',5,...
% % %         'color',[0.4 0.4 0.8],'markerfacecolor',[0.4 0.4 0.8])
% % %     hold on
% % %     set(gca,'xlim',[-0.5 1],'ylim',[-0.1 0.2])
% % %     h = lsline;
% % %     set(h,'color',[0.3 0.3 0.9],'linewidth',1.5)
% % %     plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--')
% % %     plot([0 0],get(gca,'ylim'),'color','k','linestyle','--')
% % %     xlabel('Hallway remapping')
% % %     ylabel('Room remapping')
% % %     title('Saline')
% % %     [r pval] = corr(tmp2.tmp{1},tmp{1});
% % %     text(0.4,0.16,sprintf('r = %0.3f\np = %0.3f',r,pval),...
% % %         'fontname','arial','fontsize',9,'fontweight','bold');
% % %     set(gca,'fontname','arial','fontsize',9,'fontweight','bold')  
% % %     axis square
% % %     subplot(1,2,2)
% % %     plot(tmp2.tmp{2},tmp{2},'linestyle','none','marker','o','markersize',5,...
% % %         'color',[0.4 0.4 0.8],'markerfacecolor',[0.4 0.4 0.8])
% % %     hold on
% % %     h = lsline;
% % %     set(h,'color',[0.3 0.3 0.9],'linewidth',1.5)
% % % %     set(gca,'xlim',[-0.5 1],'ylim',[-0.1 0.2])
% % %     plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--')
% % %     plot([0 0],get(gca,'ylim'),'color','k','linestyle','--')
% % %     xlabel('Hallway remapping')
% % %     ylabel('Room remapping')
% % %     title('CNO')
% % %     [r pval] = corr(tmp2.tmp{2},tmp{2});
% % %     text(0.4,-15,sprintf('r = %0.3f\np = %0.3f',r,pval),...
% % %         'fontname','arial','fontsize',9,'fontweight','bold');
% % %     set(gca,'fontname','arial','fontsize',9,'fontweight','bold')
% % %     axis square
    
    
    root = 'Plots/Summary/TwoSmall_CA3/';
    root = 'Plots/Summary/TwoSmall_CA1/';
    root = 'Plots/Summary/DREADDs_TwoSmall/';
    
    figure(2)
    saveFig(gcf,[root 'PVRemapping'],[{'tiff'} {'pdf'}]);
    figure(4)
    saveFig(gcf,[root 'PVRemapping_Animalwise'],[{'tiff'} {'pdf'}]);
    figure(6)
    saveFig(gcf,[root 'HallwayVSRoomRemapping'],[{'tiff'} {'pdf'}]);
    
%     sVals = repmat({[]},[2 nanmax(pairNum(:))]);
%     for group = 1:2
%         for i = 1:nanmax(pairNum(:))
%             sVals{group,i} = tmp{group}(pairNum(~isnan(pairNum(:,group)),group)==i);
%         end
%     end
    
%     figure(3)
%     mkGraph(sVals)
    
    
%     figure(3)
%     scatter(pairNum(~isnan(pairNum(:,1)),1),tmp{1})
%     lsline
%     set(gca,'xlim',[0 8])
%     hold on
%     scatter(pairNum(~isnan(pairNum(:,2)),2),tmp{2})
%     lsline
%     plot([0 8],[0 0],'color','k','linestyle','--')
%     xlabel('Trial pair number')
%     ylabel('Different - same similarity')
%     set(gca,'ylim',[-0.25 0.05])
    
%     [h p ci tstat] = ttest(tmp{1},tmp{2})
%     root = 'Plots/Summary/DREADDs_TwoSmall/';
% %     figure(1)
% %     saveFig(gcf,[root 'RDM'],'pdf');
%     figure(2)
%     saveFig(gcf,[root 'Similarity'],'pdf');
%     
%     outP = ['Stats_TwoBig.txt'];
%     fid = fopen(outP,'w');
%     fprintf(fid,'\t\t\tSIMILARITY\n');
%     for i = 1:2
%         for j = i+1:2
%             [h p ci tstat] = ttest(v(:,i),v(:,j));
%             fprintf(fid,['\n' num2str(i) ' to ' num2str(j) ':  ']);
%             fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
%         end
%     end
%     
%     fclose all
end

function [within across] = help_getfieldshifts(map)
    locDiff = nan(length(map(1,1,:,1)),8);

    [xg yg] = meshgrid(1:length(map(1,:,1,1)),1:length(map(:,1,1,1)));
    for k = 1:length(map(1,1,:,1))
        x = repmat({[]},[1 4]);
        y = repmat({[]},[1 4]);
        for q = 1:4
            %%% Compute based on max
            [x{q} y{q}] = find(map(:,:,k,q)== ...
                repmat(nanmax(nanmax(map(:,:,k,q),[],1),[],2),[size(map(:,:,1,1))]));
            %%% Compute based on COM
%             tmp = map(:,:,k,q);
%             tmp = tmp./nansum(tmp(:));
%             x{q} = nansum(tmp(:).*xg(:));
%             y{q} = nansum(tmp(:).*yg(:));
        end
        locDiff(k,:) = [cellfun(@nanmedian,x) cellfun(@nanmedian,y)];
    end
    locDiff = locDiff.*2.5; %convert to cm

    within = nanmean([sqrt([locDiff(:,1)-locDiff(:,3)].^2 + [locDiff(:,5)-locDiff(:,7)].^2) ...
        sqrt([locDiff(:,2)-locDiff(:,4)].^2 + [locDiff(:,6)-locDiff(:,8)].^2)],2);
    across = nanmean([sqrt([locDiff(:,1)-locDiff(:,4)].^2 + [locDiff(:,5)-locDiff(:,8)].^2) ...
        sqrt([locDiff(:,2)-locDiff(:,3)].^2 + [locDiff(:,6)-locDiff(:,7)].^2)],2);
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end









































