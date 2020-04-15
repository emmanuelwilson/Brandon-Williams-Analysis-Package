function locXRemap(paths)
    warning off all
    
    crossAllComp = [ {[]} {[]} {[]} {[]}];
    abm = [{[]} {[]}];
    
    labels = [{'Saline'} {'CNO'}];
    allIValDiffs = repmat({[]},[length(paths) 1]);
    avi = [{[]} {[]}];
    rr = [{[]} {[]}];
    doRot = [0:1:90];
    binSizeX = 4;
    binSizeY = 6;
    bbm = repmat({[]},[binSizeX binSizeY 2]);
    for p = paths'
        group = find(ismember(labels,p{1}(find(ismember(p{1},'_'),1,'last')+1:end-4)));
        if isempty(group)
            group = 1;
        end
%         if group~=2
%             continue
%         end
        
        s = load(p{1});
        
        
        isGood = s.processed.splithalf.roomXdoors.p<=0.05;
        gT = s.processed.trace(isGood,:);
        [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        
        fprintf(['\n\t' p{1} ':  ' num2str(nansum(isGood))])
        cent = s.calcium.Centroids(isGood,:);
        a = nanmean(gT(:,isInRoom&isMostRecent(1,:)),2);
        b = nanmean(gT(:,isInRoom&isMostRecent(2,:)),2);
        ival = [abs(a-b)./nanmax(a,b)];
        
        ival = nanmean(gT,2)./nanstd(gT,[],2); %tSNR
        
        half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;        
        
%         ival = s.processed.similarity.roomXdoor.ival;
        crossAllComp{(group-1).*2+1} = [crossAllComp{(group-1).*2+1}; ...
            corr(cent(:,1),ival)];
        crossAllComp{(group-1).*2+2} = [crossAllComp{(group-1).*2+2}; ...
            corr(cent(:,2),ival)];
        
%         val = nan(1,length(doRot));
%         for i = doRot
%             tc = [[cosd(i) -sind(i); sind(i) cosd(i)]*cent']';
%             val(i==doRot) = corr(tc(:,2),ival);
%         end
%         
%         rr{group} = [rr{group}; val];
        
        binMap = nan([binSizeX binSizeY]);
        cent = bsxfun(@minus,cent,nanmin(cent,[],1));
        cent = bsxfun(@rdivide,cent,nanmax(cent,[],1)+1);
        avi{group} = [avi{group}; cent ival];
        cent = floor([cent(:,1).*binSizeX cent(:,2).*binSizeY])+1;
        for ui = unique(cent,'rows')'
            binMap(ui(1),ui(2)) = nanmean(ival(ismember(cent,ui','rows'),:));
            bbm{ui(1),ui(2),group} = [bbm{ui(1),ui(2),group}; ival(ismember(cent,ui','rows'),:)];
        end
        abm{group} = cat(3,abm{group},binMap);
    end
    close all
    figure(2)
    mkGraph(crossAllComp(~cellfun(@isempty,crossAllComp)));
    drawnow
    
    figure(1)
    set(gcf,'position',[50 50 600 800])
    subplot(2,1,1)
    imagesc(nanmean(abm{1},3))
    set(gca,'ydir','normal')
    colormap hot
    axis off
    colorbar
    axis equal
    title('Saline')
    caxis([0.2 0.4])
    subplot(2,1,2)
    imagesc(nanmean(abm{2},3))
    set(gca,'ydir','normal')
    colormap hot
    axis off
    colorbar
    axis equal
    title('CNO')
    caxis([0.2 0.4])
    
%     pvals = nan(binSizeX,binSizeY);
%     for i = 1:binSizeX
%         for j = 1:binSizeY
%             [h pval ci tstat] = ttest2(bbm{i,j,1},bbm{i,j,2});
%             pvals(i,j) = pval;
%         end
%     end
    
%     [h pval ci tstat] = ttest(permute(abm{1},[3 2 1]),permute(abm{2},[3 2 1]));
%     pval = permute(pval,[3 2 1]);
%     
    root = 'Plots/Summary/DREADDs_TwoSmall/';
    figure(1)
    saveFig(gcf,[root 'RemappingXAnatomy_Heatmap'],[{'tiff'} {'pdf'}]);
    figure(2)
    saveFig(gcf,[root 'RemappingXAnatomy_Correlations'],[{'tiff'} {'pdf'}]);
end















