function plotMaps(paths)
    for p = paths'
        s = load(p{1});
        
        m = mkHDMaps(s.pos.hd,s.processed.trace,[]);
        doK = [5 5];
        for part = 0:floor(length(m)/prod(doK))

            figure(1)
            set(gcf,'position',[50 50 1000 1000])
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(m(1,:))
                    break
                end
                subplot(doK(1),doK(2),k)
                polar([0:(2.*pi./(length(m(:,1))-1)):2.*pi 0]',...
                    m([1:end 1],part.*prod(doK)+k))
                axis equal
                axis off
            end
            
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/HDMaps/Undifferentiated/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'tiff')
            saveFig(gcf,outP,'pdf')
            close all
            drawnow
        end
        
        [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        
        m1 = mkHDMaps(s.pos.hd,s.processed.trace,isMostRecent(1,:)&isInRoom);
        m2 = mkHDMaps(s.pos.hd,s.processed.trace,isMostRecent(2,:)&isInRoom);
        doK = [5 5];
        for part = 0:floor(length(m)/prod(doK))

            figure(1)
            set(gcf,'position',[50 50 1000 1000])
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(m(1,:))
                    break
                end
                subplot(doK(1),doK(2),k)
                polar([0:(2.*pi./(length(m(:,1))-1)):2.*pi 0]',...
                    m1([1:end 1],part.*prod(doK)+k));
                hold on
                polar([0:(2.*pi./(length(m(:,1))-1)):2.*pi 0]',...
                    m2([1:end 1],part.*prod(doK)+k));
                axis equal
                axis off
            end
            
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/HDMaps/Differentiated/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'tiff')
            saveFig(gcf,outP,'pdf')
            close all
            drawnow
        end
    end
end