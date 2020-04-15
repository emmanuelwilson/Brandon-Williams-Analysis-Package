function mapsXDoors5_parts(paths)
    crossAllComp = [];
    crossDataAmount = [];
    labels = [];
    groups = [{'Saline'} {'CNO'}];
    for pi = length(paths):-1:1'
        p = paths(pi);
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        
%         gT = s.processed.trace(s.processed.isAligned,:);
        gT = s.processed.trace(s.processed.splithalf.p<=0.05,:);
%         gT = s.processed.trace(s.processed.splithalf.within.p<=0.05,:);
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
        doRange = [0:4:24];
%         doRange = [0:7.5:15];
%         doK = [5 9 10 13 31 34]; % 5 9 10 13 31 34
        rangeComp = nan(2,length(doRange));
        dataAmount = nan(1,length(doRange));
        av = nan(2.*(length(doRange)-1));
        for curRI = 1:length(doRange)-1
            indexSinceIn = oIndexSinceIn;
            indexSinceIn = nanmin(indexSinceIn).*(1./30) >= doRange(curRI) & ...
                nanmin(indexSinceIn).*(1./30) < doRange(curRI+1);
            dataAmount(curRI) = nansum(indexSinceIn);
            if ~any(indexSinceIn)
                continue
            end
            for curRJ = 1:length(doRange)-1
                indexSinceIn2 = oIndexSinceIn;
                indexSinceIn2 = nanmin(indexSinceIn2).*(1./30) >= doRange(curRJ) & ...
                    nanmin(indexSinceIn2).*(1./30) < doRange(curRJ+1);
                dataAmount(curRI) = nansum(indexSinceIn);
                if ~any(indexSinceIn)
                    continue
                end
                for doorA = 1:2
                    for doorB = 1:2
                        if (doorA-1).*(length(doRange)-1)+curRI <= (doorB-1).*(length(doRange)-1)+curRJ
                            continue
                        end
                        [map samp val] = getMatchedSamplingMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),...
                            [isMostRecent(doorA,isInRoom) & indexSinceIn(1,isInRoom)],...
                            [isMostRecent(doorB,isInRoom) & indexSinceIn2(1,isInRoom)]);
                        av((doorA-1).*(length(doRange)-1)+curRI,...
                            (doorB-1).*(length(doRange)-1)+curRJ) = val;
                    end
                end
            end
        end
        tav = av;
        tav = nanmean(cat(3,tav,tav'),3);
        tav(1:(length(doRange)-1).*2+1:end) = 1;
        pts = mdscale(1-tav,2);
        for k = 1:(length(doRange)-1)
            hold on
            plot([pts(k,1) pts((length(doRange)-1)+k,1)],...
                [pts(k,2) pts((length(doRange)-1)+k,2)],'color','k',...
                'linewidth',1)
            scatter(pts(k,1),pts(k,2),k.*15,'r','filled')
            scatter(pts((length(doRange)-1)+k,1),...
                pts((length(doRange)-1)+k,2),k.*15,'b','filled')
        end
    end
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end









































