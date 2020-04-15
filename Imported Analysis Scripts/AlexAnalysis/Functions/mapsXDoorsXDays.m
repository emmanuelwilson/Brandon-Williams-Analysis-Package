function mapsXDoorsXDays(folder)

    mice = dir(folder);
    mice = {mice(3:end).name};
    
    similarity = repmat({[]},[1 length(mice)]);
    for m = mice
        mp = [folder '/' m{1}];
        sessions = dir(mp);
        sessions = {sessions(3:end).name};
        
        % order sessions by date
        clear dates
        for j = 1:length(sessions)
            dates(j) = datetime(str2num(sessions{j}(1:2)),str2num(sessions{j}(4:5)),str2num(sessions{j}(7:8)));
        end
        [a b] = sort(dates);
        sessions = sessions(b);
        
        com = repmat({[]},[1 2]);
        uni = repmat({[]},[1 2]);
        isInRoom = repmat({[]},[1 2]);
        
        for dci = 1:length(sessions)
            s1 = load([mp '/' sessions{dci}]);
            
            [isIn isMostRecent] = isInROI(s1.processed.p,s1.processed.roi.door);
        
            com{1} = isMostRecent(1,:)|isMostRecent(3,:)|isMostRecent(5,:);
            uni{1} = isMostRecent(2,:)|isMostRecent(4,:)|isMostRecent(6,:);

            isInRoom{1} = isInROI(s1.processed.p,s1.processed.roi.room);
            
            
            for dcj = dci+1:length(sessions)
                
                s2 = load([mp '/' sessions{dcj}]);
                
                [isIn isMostRecent] = isInROI(s2.processed.p,s2.processed.roi.door);
        
                com{2} = isMostRecent(1,:)|isMostRecent(3,:)|isMostRecent(5,:);
                uni{2} = isMostRecent(2,:)|isMostRecent(4,:)|isMostRecent(6,:);

                isInRoom{2} = isInROI(s2.processed.p,s2.processed.roi.room);


                doSize = zeros(1,3);
                for room = 1:3
                    m1 = mkTraceMaps(s1.processed.p(:,[com{1} & isInRoom{1}(room,:)]),...
                        s1.processed.trace(s1.processed.isAligned,[com{1} & isInRoom{1}(room,:)]),[]);
                    m2 = mkTraceMaps(s1.processed.p(:,[uni{1} & isInRoom{1}(room,:)]),...
                        s1.processed.trace(s1.processed.isAligned,[uni{1} & isInRoom{1}(room,:)]),[]);
                    m3 = mkTraceMaps(s2.processed.p(:,[com{2} & isInRoom{2}(room,:)]),...
                        s2.processed.trace(s2.processed.isAligned,[com{2} & isInRoom{2}(room,:)]),[]);
                    m4 = mkTraceMaps(s2.processed.p(:,[uni{2} & isInRoom{2}(room,:)]),...
                        s2.processed.trace(s2.processed.isAligned,[uni{2} & isInRoom{2}(room,:)]),[]);
                    doSize = nanmax([doSize; size(m1); size(m2); size(m3); size(m4)]);
                end

                comMaps = nan([(doSize) 3 2]);
                uniMaps = nan([(doSize) 3 2]);

                for room = 1:3
                    comMaps(:,:,:,room,1) = mkTraceMaps(s1.processed.p(:,[com{1} & isInRoom{1}(room,:)]),...
                        s1.processed.trace(s1.processed.isAligned,[com{1} & isInRoom{1}(room,:)]),[],doSize(1:2));
                    uniMaps(:,:,:,room,1) = mkTraceMaps(s1.processed.p(:,[uni{1} & isInRoom{1}(room,:)]),...
                        s1.processed.trace(s1.processed.isAligned,[uni{1} & isInRoom{1}(room,:)]),[],doSize(1:2));
                    comMaps(:,:,:,room,2) = mkTraceMaps(s2.processed.p(:,[com{2} & isInRoom{2}(room,:)]),...
                        s2.processed.trace(s2.processed.isAligned,[com{2} & isInRoom{2}(room,:)]),[],doSize(1:2));
                    uniMaps(:,:,:,room,2) = mkTraceMaps(s2.processed.p(:,[uni{2} & isInRoom{2}(room,:)]),...
                        s2.processed.trace(s2.processed.isAligned,[uni{2} & isInRoom{2}(room,:)]),[],doSize(1:2));
                end
                
                

                comMaps = cat(4,comMaps,uniMaps);
                
                allComp = nan(6,6);
                for ri = 1:6
                    for rj = 1:6
                        m1 = comMaps(:,:,:,ri,1);
                        m2 = comMaps(:,:,:,rj,2);
                        isGood = ~(isnan(m1)|isnan(m2));
                        allComp(ri,rj) = corr(m1(isGood),m2(isGood));
                    end
                end
                
                similarity{ismember(mice,m)} = cat(3,...
                    similarity{ismember(mice,m)},allComp);                
            end
        end
    end 
    
    imagesc(nanmean(cat(3,similarity{:}),3))
    colormap jet
    colorbar
    caxis([-0.1 0.35])
    axis equal
end



































































