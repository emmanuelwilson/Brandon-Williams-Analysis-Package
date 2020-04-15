function mapsXDoors5_wholemaps(paths)
    crossAllComp = [{[]} {[]}];
    labels = [{'Saline'} {'CNO'}];
    allIValDiffs = repmat({[]},[length(paths) 1]);
    for p = paths'
        s = load(p{1});
        group = find(ismember(labels,p{1}(find(ismember(p{1},'_'),1,'last')+1:end-4)));
        if isempty(group)
            group = 1;
        end
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
        gT = s.processed.trace(s.processed.splithalf.roomXdoors.p<=0.05,:);
%         gT = s.processed.trace(s.processed.splithalf.within.p<=0.05,:);
%         gT = s.processed.trace;

        fprintf(['\n\t' p{1} ':  ' num2str(length(gT(:,1)))])

        [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        a = nanmin(indexSinceIn).*(1./30);
        thresh = -1;

        doSize = zeros(1,3);
        for room = 1:1
            m1 = mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]));
            doSize = nanmax([doSize; size(m1)]);
        end
        
        allMasks = repmat({[]},[1 2]);
        for i = 1:2
            allMasks{i} = [isMostRecent(i,isInRoom)];
        end
        
        [map samp allComp ivals] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks(1:2));
        
        
%         crossRateComp = cat(3,crossRateComp,rateComp);
        crossAllComp{group} = cat(2,crossAllComp{group},allComp(2));
%         crossAllComp{group} = cat(3,crossAllComp{group},nanmean(ivals(1:2,3:4,:),3));
    end

    
    %% Average within animal
    piece = [];
    ag = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        ag = [ag; find(ismember(labels,paths{i}(find(ismember(paths{i},'_'),1,'last')+1:end-4)))];
    end
    upiece = unique(piece);
    pairNum = nan(length(piece),2);
    count = zeros(1,length(upiece));
    for i = 1:length(piece)
        count(ismember(upiece,piece(i))) = count(ismember(upiece,piece(i)))+1;
        pairNum(i,ag(i)) = ceil(count(ismember(upiece,piece(i)))./2);
    end

    
%     tmp = cellfun(@nanmedian,allIValDiffs);
%     tmp(isnan(pairNum(:,1))) = [];
%     mkGraph([{tmp(1:4)} {tmp(5:10)} {tmp(11:16)}]);
    
    %%%%%%%%%
    
    av = crossAllComp;
    tmp = av;
%     tmp = cellfun(@diff,av,[{[]} {[]}],[{2} {2}],'uniformoutput',false);
%     tmp = cellfun(@rdivide,tmp,[{av{1}(:,1)} {av{2}(:,1)}],'uniformoutput',false);
%     tmp2 = tmp{1}-tmp{2};
%     mkGraph(tmp)
    
    figure(4)
    animalNormed = repmat({[]},[length(upiece) 2]);
    for i = 1:length(upiece)
        for j = 1:2
            animalNormed{i,j} = tmp{j}(ismember(piece(ag==j),upiece(i)))';
        end
    end
%     hist(tmp2,[-0.15:0.015:0.15])
    mkGraph(animalNormed');
    
    [h p ci tstat] = ttest(animalNormed{3,1},animalNormed{3,2})
    
    sVals = repmat({[]},[2 nanmax(pairNum(:))]);
    for group = 1:2
        for i = 1:nanmax(pairNum(:))
            sVals{group,i} = tmp{group}(pairNum(~isnan(pairNum(:,group)),group)==i);
        end
    end
    
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









































