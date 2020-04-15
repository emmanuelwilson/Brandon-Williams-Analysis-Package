function alignPairwiseSessions(folder)
    mice = dir(folder);
    mice = {mice(3:end).name};
    for m = mice
        mp = [folder '/' m{1}];
        sessions = dir(mp);
        sessions = {sessions(3:end).name};
        
        if isempty(sessions)
            continue
        end
        
        % order sessions by date
        clear dates
        for j = 1:length(sessions)
            dates(j) = datetime(str2num(sessions{j}(1:2)),str2num(sessions{j}(4:5)),str2num(sessions{j}(7:8)));
        end
        [a b] = sort(dates);
        sessions = sessions(b);
        
        
        alignmentMap = repmat({[]},[length(sessions) length(sessions)]);
        for j = 1:2:length(sessions)
            for k = j+1
            
                if j+1 > length(sessions) %%% In case odd number of sessions
                    continue
                end

                for si = [j k]
                    ref = load([mp '/' sessions{si}]);                
%                     prepped = msExtractSFPs(ref.calcium);
                    prepped = ref.calcium.SFPs .* ...
                        bsxfun(@gt,ref.calcium.SFPs,0.5.*nanmax(nanmax(ref.calcium.SFPs,[],1),[],2));

                    outP = ['SegmentsForAlignment/' sessions{si}];
                    checkP(outP)
                    save(outP,'prepped'); 
                end

                map = registerCells(['SegmentsForAlignment']);
                rmdir(['SegmentsForAlignment'],'s');
                close all
                drawnow

                isReg = map(all(map~=0,2),:);
                alignmentMap{j,k} = isReg;
            end
        end
        ref = load([mp '/' sessions{1}]);
        ref.alignmentMap = alignmentMap;
        save([mp '/' sessions{1}],'-struct','ref','-v7.3');

    end
end