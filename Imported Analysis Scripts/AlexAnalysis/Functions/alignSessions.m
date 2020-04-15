function alignSessions(folder)
    mice = dir(folder);
    mice = {mice(3:end).name};
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
        
        for j = 1:1:length(sessions)
            ref = load([mp '/' sessions{j}]);
            
            % prep mat
            prepped = ref.calcium.segments.*repmat(ref.calcium.frameMax,[1 1 length(ref.calcium.trace(1,:))]);
            prepped = double(prepped) ./ repmat(nansum(nansum(prepped,1),2),[size(prepped(:,:,1)) 1]);
            prepped = permute(prepped,[3 1 2]);
            
            outP = [mp '/SegmentsForAlignment/' sessions{j}];
            checkP(outP)
            save(outP,'prepped');
        end
        map = registerCells([mp '/SegmentsForAlignment']);
        rmdir([mp '/SegmentsForAlignment'],'s');
        close all
        drawnow
        
        %%% Save successfully registered cells
        isReg = map(all(map~=0,2),:);
        for j = 1:length(sessions)
            ref = load([mp '/' sessions{j}]);
            
            isNotReg = ~ismember([1:length(ref.processed.trace(:,1))]',isReg(:,j));
            ref.processed.trace = [ref.processed.trace(isReg(:,j),:); ref.processed.trace(isNotReg,:)];
            ref.processed.isAligned = [true(length(isReg(:,1)),1); false(sum(isNotReg),1)]; 
            ref.processed.alignmentMap = map;
            
            save([mp '/' sessions{j}],'-struct','ref','-v7.3');
        end
    end
end