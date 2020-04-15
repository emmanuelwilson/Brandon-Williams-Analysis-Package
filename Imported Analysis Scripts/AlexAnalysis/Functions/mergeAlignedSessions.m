function dat = mergeAlignedSessions(folder,outFolder)
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
        
        s = load([mp '/' sessions{1}]);
        registration = s.alignmentMap;
        
        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                if ~isempty(registration{si,sj})
                    s1 = load([mp '/' sessions{si}]);
                    s2 = load([mp '/' sessions{sj}]);
                    newS = s1;
                    newS.processed.p = [s1.processed.p s2.processed.p];
                    newS.processed.trace = [s1.processed.trace(registration{si,sj}(:,1),:) ...
                        s2.processed.trace(registration{si,sj}(:,2),:)];
                    newS.processed.roi = [s1.processed.roi s2.processed.roi];
                    
                    outP = [outFolder '/MERGED_' sessions{si}(1:end-4) ...
                        '_' sessions{sj}(1:end-4) '.mat'];
                    checkP(outP);
                    save(outP,'-struct','newS','-v7.3');
                end
            end
        end        
    end    
end