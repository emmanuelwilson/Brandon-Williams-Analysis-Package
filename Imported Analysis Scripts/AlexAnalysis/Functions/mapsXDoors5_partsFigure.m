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
        doRange = [0:3:18];
%         doRange = [0:7.5:15];
        doK = [5 9 10 13 31 34]; % 5 9 10 13 31 34
        rangeComp = nan(2,length(doRange));
        dataAmount = nan(1,length(doRange));
        figure
        set(gcf,'position',[50 50 250.*length(doRange)-1 length(doK).*175])
        for curR = 1:length(doRange)-1
            indexSinceIn = oIndexSinceIn;
%             indexSinceIn = nanmin(indexSinceIn).*(1./30) >= doRange(curR) & ...
%                 nanmin(indexSinceIn).*(1./30) < doRange(curR+1);
            indexSinceIn = nanmin(indexSinceIn).*(1./30) >= doRange(curR);
            dataAmount(curR) = nansum(indexSinceIn);
            if ~any(indexSinceIn)
                continue
            end
            m1 = mkTraceMaps(s.processed.p(:,isInRoom),...
                gT(:,isInRoom),isMostRecent(1,isInRoom)&indexSinceIn(isInRoom));
            m2 = mkTraceMaps(s.processed.p(:,isInRoom),...
                gT(:,isInRoom),isMostRecent(2,isInRoom)&indexSinceIn(isInRoom));
            
            m = cat(2,m1,nan(size(m1(:,1,:))),m2);
            
            for k = 1:length(doK)
                subplot(length(doK),length(doRange)-1,(k-1).*(length(doRange)-1)+curR)
                imagesc(m(:,:,doK(k)))
                set(gca,'ydir','normal')
                colormap jet
                alpha(double(~isnan(m(:,:,doK(k)))))
                axis equal
                axis off   
            end
            
          
        end
        slashInds = find(ismember(p{1},'/'));
        outP = ['Plots/DifferentiatedMaps_ByParts/' p{1}(slashInds+1:end-4)];
        saveFig(gcf,outP,'tiff')
        saveFig(gcf,outP,'pdf')
        close all
        drawnow
    end
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end









































