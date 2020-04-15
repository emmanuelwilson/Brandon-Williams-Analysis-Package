function alignAllSessions(paths)
    clc
    fprintf('\n')
    %%% Reliability constaint
    
    %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
    
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}) '\n'])
        
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);      
        
        alignmentMap = repmat({[]},[length(sessions) length(sessions)]);

        for si = 1:length(sessions)
            ref = load(sessions{si});                
            prepped = msExtractSFPs(ref.calcium);
            outP = ['SegmentsForAlignment/' num2str(si)];
            checkP(outP)
            save(outP,'prepped'); 
        end

        try
            map = registerCells(['SegmentsForAlignment']);
        catch
            map = [];
        end
        rmdir(['SegmentsForAlignment'],'s');
        close all
        drawnow

        isReg = map(all(map~=0,2),:);
        alignmentMap = isReg;

        ref = load(sessions{1});
        ref.alignment.allAlignmentMap = alignmentMap;
        save(sessions{1},'-struct','ref','-v7.3');

    end
end