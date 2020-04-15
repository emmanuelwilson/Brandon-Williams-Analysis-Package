function dat = getTrialData(folder)
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
            gT = s.processed.trace;
            
            if isfield(s.processed,'splithalf')
                dat(ismember(mice,m)).shr{si} = s.processed.splithalf.p;
            end
            dat(ismember(mice,m)).trace{si} = gT;
            dat(ismember(mice,m)).path{si} = s.processed.p;
            

            dat(ismember(mice,m)).maps.overall{si} = mkTraceMaps(s.processed.p,gT,[],[11 21]);
        end        
    end

    for mi = length(dat):-1:1
        if isempty(dat(mi).maps.overall)
            dat(mi) = [];
        end
    end

    newDat.dat = dat;
    save('combinedData_OpenField','-struct','newDat','-v7.3');
    
end