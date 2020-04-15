function plotRoomMaps(paths)

    doColors = [1 1 1; ...
        0.35 0.0 0.9; ...
        0.0 0.9 0.9];
    
    cstack = [ ... %repmat(permute(doColors(1,:),[3 1 2]),[12 36]); ...
        repmat(permute(doColors(2,:),[3 1 2]),[12 36]); ...
        repmat(permute(doColors(3,:),[3 1 2]),[12 36])];
    
    for p = paths'
        s = load(p{1});
        gT = s.processed.trace;
        try
            [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
        catch
            continue
        end
        com = isMostRecent(1,:)|isMostRecent(3,:)|isMostRecent(5,:);
        uni = isMostRecent(2,:)|isMostRecent(4,:)|isMostRecent(6,:);
        [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
        clear maps
        for room = 1:3
            maps.rooms{room} = mkTraceMaps( ...
                s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),[],[12 12]);

            maps.common{room} = ...
                mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),com(isInRoom(room,:)),[12 12]);

            maps.unique{room} = ...
                mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                gT(:,[isInRoom(room,:)]),uni(isInRoom(room,:)),[12 12]);
        end
        
        m = nan(24,36,length(maps.rooms{1}(1,1,:)));
        for r = 1:3
%             m(1:12,(r-1).*12+1:(r).*12,:) = maps.rooms{r};
            m([13:24]-12,(r-1).*12+1:(r).*12,:) = maps.common{r};
            m([25:36]-12,(r-1).*12+1:(r).*12,:) = maps.unique{r};
        end
        
        doK = [12 12];
        for part = 0:floor(length(m)/prod(doK))

            figure(1)
            set(gcf,'position',[50 50 1400 900])
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(m(1,1,:))
                    break
                end
                subplot(doK(1),doK(2),k)
                tm = m(:,:,part.*prod(doK)+k);
                tm = tm./nanmax(tm(:));
%                 image(bsxfun(@times,tm,cstack));
                imagesc(m(:,:,part.*prod(doK)+k))
                colormap('jet')
                alpha(double(~isnan(m(:,:,part.*prod(doK)+k))))
%                 title(num2str(nanmax(nanmax(m(:,:,part.*prod(doK)+k)))))
%                 drawnow
%                 set(gca,'ydir','normal')
                axis equal
                axis off
            end
            
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/TraceFieldRoomMaps_Differentiated/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'tiff')
            saveFig(gcf,outP,'pdf')
            close all
            drawnow
        end
    end
end