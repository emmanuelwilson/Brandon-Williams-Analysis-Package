function dat = getTrialData(folder)
    mice = dir(folder);
    mice = {mice(3:end).name};
    
%     mice = mice(4);
    
    maps = struct('overall',[],'rooms',[],'common',[],'unique',[]);
    entrances = struct('common',[],'unique',[]);
    dat = struct('maps',maps,'shr',[],'trace',[],'path',[],'isAligned',[], ...
        'entrances',entrances,'mostrecent',entrances);
    dat = repmat(dat,[1 length(mice)]);
    for m = mice
        mp = [folder '/' m{1}];
        sessions = dir(mp);
        sessions = {sessions(3:end).name};
        
        % order sessions by date
%         clear dates
        dates = datetime([], 'ConvertFrom', 'datenum');
        for j = 1:length(sessions)
            dates(j) = datetime(str2num(sessions{j}(1:2)),str2num(sessions{j}(4:5)),str2num(sessions{j}(7:8)));
        end
        [a b] = sort(dates);
        sessions = sessions(b);
        
        for si = 1:length(sessions)
            s = load([mp '/' sessions{si}]);

            if isfield(s,'alignmentMap')
                dat(ismember(mice,m)).registration = s.alignmentMap;
            end
            
%             gT = s.processed.trace(s.processed.splithalf.p<=0.05,:);
            gT = s.processed.trace;
            
            if isfield(s.processed,'splithalf')
                dat(ismember(mice,m)).shr{si} = s.processed.splithalf.p;
            end
            
            dat(ismember(mice,m)).trace{si} = gT;
            dat(ismember(mice,m)).path{si} = s.processed.p;
            
            [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
        
            com = isMostRecent(1,:)|isMostRecent(3,:)|isMostRecent(5,:);
            uni = isMostRecent(2,:)|isMostRecent(4,:)|isMostRecent(6,:);
            
            [isInRoom a b] = isInROI(s.processed.p,s.processed.roi.room);
            
            dat(ismember(mice,m)).mostrecent.common{si} = com;
            dat(ismember(mice,m)).mostrecent.unique{si} = uni;
            dat(ismember(mice,m)).isinroom{si} = isInRoom;
%             com = com & nanmin(indexSinceIn([1 3 5],:),[],1) > 30; % toss the first 30 frames (1 sec) when entering a room
%             uni = uni & nanmin(indexSinceIn([2 4 6],:),[],1) > 30; % toss the first 30 frames (1 sec) when entering a room

            [isInRoom a b] = isInROI(s.processed.p,s.processed.roi.room);
            
            %%% If just exited the room and came back in, ignore the
            %%% exit/entry????????????
            
            
            isRoomEntry = diff(any(isInRoom),[],2)>0;
            dat(ismember(mice,m)).entrances.common(si) = nanmean(com(isRoomEntry));
            dat(ismember(mice,m)).entrances.unique(si) = nanmean(uni(isRoomEntry));
            
            dat(ismember(mice,m)).maps.overall{si} = mkTraceMaps(s.processed.p,gT,[],[20 40]);
            
            for room = 1:3
                
                dat(ismember(mice,m)).maps.rooms{room,si} = mkTraceMaps( ...
                    s.processed.p(:,[isInRoom(room,:)]),...
                    gT(:,[isInRoom(room,:)]),[],[12 12]);
                
                dat(ismember(mice,m)).maps.common{room,si} = ...
                    mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                    gT(:,[isInRoom(room,:)]),com(isInRoom(room,:)),[12 12]);
                
                dat(ismember(mice,m)).maps.unique{room,si} = ...
                    mkTraceMaps(s.processed.p(:,[isInRoom(room,:)]),...
                    gT(:,[isInRoom(room,:)]),uni(isInRoom(room,:)),[12 12]);
            end
        end        
    end

    for mi = length(dat):-1:1
        if isempty(dat(mi).maps.overall)
            dat(mi) = [];
        end
    end

    newDat.dat = dat;
    save('combinedData_Exp1','-struct','newDat','-v7.3');
end