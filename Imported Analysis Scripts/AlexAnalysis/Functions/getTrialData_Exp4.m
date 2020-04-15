function dat = getTrialData_Exp(folder)
    mice = dir(folder);
    mice = {mice(3:end).name};
    
%     mice = mice(4);
    
    maps = struct('overall',[],'rooms',[],'common',[],'unique',[]);
    dat = struct('maps',maps,'shr',[],'trace',[],'path',[],'isAligned',[],'registration',[]);
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
            
            %%% Eliminate paths that were not initiated at the start of the
            %%% long path (BELOW)
            
            [isIn isMostRecent] = isInROI(s.processed.p,s.processed.roi.door);
            
%             inds = find(diff(isMostRecent(4,:))>0);
%             eliminate = inds(~isMostRecent(3,inds));
%             
%             for i = 1:length(eliminate)
%                 stop = find(~isMostRecent(4,eliminate(i)+1:end),1,'first');
%                 if isempty(stop)
%                     stop = length(isMostRecent(1,:))-eliminate(i);
%                 end
%                 isMostRecent(4,eliminate(i):eliminate(i)+stop) = false;
%             end
            isMostRecent(3,:) = isMostRecent(4,:);
            isMostRecent(4,:) = [];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isfield(s.processed,'splithalf')
                dat(ismember(mice,m)).shr{si} = s.processed.splithalf.p;
            end
%             if si ~= length(sessions)
%                 dat(ismember(mice,m)).registration(si,:,1) = s.alignmentMap(si,:);
%             end
%             dat(ismember(mice,m)).registration(1:length(s.alignmentMap(:,si)),si,2) = ...
%                 s.alignmentMap(:,si);
%             dat(ismember(mice,m)).isAligned{si} = s.processed.isAligned;
            dat(ismember(mice,m)).trace{si} = gT;
            dat(ismember(mice,m)).path{si} = s.processed.p;
            
            [isInRoom a b] = isInROI(s.processed.p,s.processed.roi.room);
            dat(ismember(mice,m)).isInRoom{si} = isInRoom;
            dat(ismember(mice,m)).door{si} = isMostRecent;
            
            [dat(ismember(mice,m)).maps.overall{si} blah] = mkTraceMaps(s.processed.p,gT,[],[24 35]);
%             if si >= 4 %%%%% SHIFT TO ALIGN BECAUSE OF THE CHANGE OF SHAPE
%                 dat(ismember(mice,m)).maps.overall{si} = ...
%                     circshift(dat(ismember(mice,m)).maps.overall{si},[0 12 0]);
%             end
            
            dat(ismember(mice,m)).maps.rooms{1,si} = mkTraceMaps( ...
                s.processed.p(:,isInRoom),...
                gT(:,isInRoom),[],[18 18]);
            
            for i = 1:3
                dat(ismember(mice,m)).maps.unique{i,si} = mkTraceMaps( ...
                    s.processed.p(:,isInRoom),...
                    gT(:,isInRoom),isMostRecent(i,isInRoom),[18 18]);
            end
        end        
    end

    for mi = length(dat):-1:1
        if isempty(dat(mi).maps.overall)
            dat(mi) = [];
        end
    end

    newDat.dat = dat;
    save('combinedData_Exp4','-struct','newDat','-v7.3');
    
end