function hallwaysXDir_5(paths)
    warning off all
    
    aiv = [{[]} {[]}];
    apfr = [{[]} {[]}];
    crossAllComp = [{[]} {[]}];
    labels = [{'Saline'} {'CNO'}];
    allIValDiffs = repmat({[]},[length(paths) 1]);
    for p = paths'
        s = load(p{1});
        group = find(ismember(labels,p{1}(find(ismember(p{1},'_'),1,'last')+1:end-4)));
        if isempty(group)
            group = 1;
        end
        
%         gT = s.processed.trace(s.processed.isAligned,:);
%         gT = s.processed.trace(s.processed.splithalf.p<=0.01,:);
%         gT = s.processed.trace(s.processed.splithalf.within.p<=0.05,:);
        gT = s.processed.trace(s.processed.splithalf.hallwayXdoors.p<=0.05,:);
%         gT = s.processed.trace;
%         gT = gT>0;

        fprintf(['\n\t' p{1} ':  ' num2str(length(gT(:,1)))])

        isIn = inpolygon(s.processed.p(1,:)',s.processed.p(2,:)',...
            s.processed.roi.hallway(:,1),s.processed.roi.hallway(:,2));        
        cp = help_collapseToLine(s.processed.p,s.processed.roi.hallway_linear);
        
        isGood = isIn; %%% Only include full movements
        isDir = false(2,length(isIn));
        oneEntrance = cp > [nanmax(cp(isIn))-nanmin(cp(isIn))]./2;
        while any(isGood)
            start = find(isGood,1,'first');
            stop = find(~isGood(start:end),1,'first')-2;
            if isempty(stop)
                stop = length(isGood)-start;
            end
%             if oneEntrance(start)~=oneEntrance(start+stop)
                if oneEntrance(start)
                    isDir(1,start:start+stop) = true;
                else
                    isDir(2,start:start+stop) = true;
                end
%             end
            isGood(start:start+stop) = false;
        end
        
        half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;

        [maps samp allComp ival] = getMatchedLinMapsNMasks(cp(isIn),gT(:,isIn),...
            [{isDir(1,isIn) & half(isIn)} {isDir(2,isIn)& half(isIn)} ...
            {isDir(1,isIn) & ~half(isIn)} {isDir(2,isIn)& ~half(isIn)}]);
        
        
        crossAllComp{group} = cat(3,crossAllComp{group},allComp(1:2,3:4));

        maps = getMatchedLinMapsNMasks(cp(isIn),gT(:,isIn),...
            [{isDir(1,isIn)} {isDir(2,isIn)}]);
        a = (permute(maps(:,:,:,1),[3 1 2]));
        b = (permute(maps(:,:,:,2),[3 1 2]));
        [blah mind] = nanmax(a,[],2);
        [blah sind] = sort(mind);
        a = a(sind,:);
        b = b(sind,:);
        normer = nanmax([a b],[],2);
        a = bsxfun(@rdivide,a,normer);
        b = bsxfun(@rdivide,b,normer);
        
        figure(2)
        colormap jet
        subplot(1,4,1)
        imagesc(a)
        axis off
        subplot(1,4,2)
        imagesc(b)
        axis off
        [blah mind] = nanmax(b,[],2);
        [blah sind] = sort(mind);
        a = a(sind,:);
        b = b(sind,:);
        subplot(1,4,3)
%         plot([a+bsxfun(@times,ones(size(a(1,:))),[1:length(a(:,1))]'.*1)]',...
%             'color',[0.25 0.0 0.0])
%         hold on
%         plot([b+bsxfun(@times,ones(size(a(1,:))),[1:length(a(:,1))]'.*1)]',...
%             'color','k')
        imagesc(a)
        axis off
        subplot(1,4,4)
        imagesc(b)
        axis off
        set(gcf,'position',[50 50 1200 800])
        axis off
        drawnow
        slashInds = find(ismember(p{1},'/'));
        outP = ['Plots/HallwayDirectionMaps/' p{1}(slashInds+1:end-4)];
        saveFig(gcf,outP,[{'tiff'} {'pdf'}])
        close all
        drawnow

        
%         aiv{group} = [aiv{group}; permute(ival(1,2,:),[3 1 2])];

        
%         close all
%         figure(3)
%         set(gcf,'position',[650 50 300 300])
%         cumHist(aiv,[-1:0.01:1])
%         drawnow
        
%         allMasks = repmat({[]},[1 4]);
%         for i = 1:2
%             allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom)];
%             allMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom)];
%         end
%         
%         map = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),{[]});
%         apfr{group} = [apfr{group}; permute(nanmax(nanmax(map,[],1),[],2),[3 2 1])];
    end
%     close all
%     figure(3)
%     set(gcf,'position',[650 50 300 300])
%     cumHist(aiv,[-1:0.01:1])

    
    close all
    
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
    
    av = [{[]} {[]}];
    for group = 1:2
    
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

        figure(2)
        set(gcf,'position',[50 450 400 200])
        subplot(1,3,group)
        mkGraph(v)
        set(gca,'ylim',[0 0.7])
        av{group} = v;
    end
	subplot(1,3,3)  
    av = cellfun(@fliplr,av,'uniformoutput',false);
    tmp = cellfun(@diff,av,[{[]} {[]}],[{2} {2}],'uniformoutput',false);
%     tmp = cellfun(@rdivide,tmp,[{av{1}(:,1)} {av{2}(:,1)}],'uniformoutput',false);
%     tmp2 = tmp{1}-tmp{2};
    mkGraph(tmp)
    save('DREADDs_TwoSmall_HallwayRemapping','tmp')
    
    
    [h p ci tstat] = ttest2(tmp{1},tmp{2})
    
    figure(4)
    set(gcf,'position',[50 350 350 250])
    animalNormed = repmat({[]},[length(upiece) 2]);
    for i = 1:length(upiece)
        for j = 1:2
            animalNormed{i,j} = tmp{j}(ismember(piece(ag==j),upiece(i)));
        end
%         if length(animalNormed{i,1})==length(animalNormed{i,2})
%             animalNormed{i,3} = animalNormed{i,1}-animalNormed{i,2};
%         end
    end
%     hist(tmp2,[-0.15:0.015:0.15])
    mkGraph(animalNormed');
    
    
    root = 'Plots/Summary/DREADDs_TwoSmall/';
    figure(2)
    saveFig(gcf,[root 'Hallway_PVRemapping'],[{'tiff'} {'pdf'}]);
    figure(4)
    saveFig(gcf,[root 'Hallway_PVRemapping_Animalwise'],[{'tiff'} {'pdf'}]);
    
    
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


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end









































