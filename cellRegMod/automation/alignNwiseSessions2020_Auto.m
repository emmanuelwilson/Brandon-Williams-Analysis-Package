%% Modified cellReg script, combination of Alex's and Emmmanuel's

function alignNwiseSessions2020_Auto(SFPsShifted,sessions,combs,path,nonRigid)
nFold = 2;
oldcd = pwd;
cd(path);
paths = getFilePaths(path,'.mat');
for i = length(paths):-1:1
    if contains(paths{i},'OriginalFiles')
        paths(i) = [];
    elseif contains(paths{i},'Segments')
        paths(i) = [];
    elseif contains(paths{i},'tempfigs')
        paths(i) = [];
    elseif contains(paths{i},'Results')
        paths(i) = [];
    end
end

iters = 10;

warning ('off','all');

clc
fprintf('\n')
%%% Reliability constaint

%% Split by animal
mkdir(path,'Results');
alignmentMap = repmat({[]},[length(combs(:,1)) 1]);
probMap = alignmentMap;
cormat = diag(max(max(combs)));    
    
%     allPrepped  = cell(length(combs(:,1)));%= repmat({[]},[1 length(sessions)]);
    for si = 1:length(combs(:,1))
        if si == 17
        end
%         fprintf(['\t\tSession:  ' sessions{combs(si,1)}(find(ismember(sessions{combs(si,1)},'/'),1,'last')+1:end) '\n'])
        ref = sessions{combs(si,1)};%,'calcium','processed');        
        preppedr = ref.SFPs .* bsxfun(@gt,ref.SFPs,0.5.*nanmax(nanmax(ref.SFPs,[],1),[],2));
        move = SFPsShifted{combs(si,1),combs(si,2)};%,'calcium','processed');        
        preppedm = move.* bsxfun(@gt,move,0.5.*nanmax(nanmax(move,[],1),[],2));
        [cellmapr,excluder] = msExtractSFPsCellReg2020(preppedr);
        [cellmapm,excludem] = msExtractSFPsCellReg2020(preppedm);
        outofFOVr.cellmap{si} = cellmapr;
        outofFOVm.cellmap{si} = cellmapm;
        outofFOVr.exclude_outofFOV{si} = excluder;
        outofFOVm.exclude_outofFOV{si} = excludem;
%         allPrepped{si,1} = permute(preppedr,[3 1 2]); %msExtractSFPs(ref.calcium);
        allPrepped = permute(preppedm,[3 1 2]); 
        if isfield(ref,'processed')
            if isfield(ref.processed,'exclude')
                excluder(find(ref.processed.exclude.SFPs) == 0) = 0;
                outofFOVr.exclude_badSFP = ref.processed.exclude.SFPs;
            end
        end
        if isfield(move,'processed')
            if isfield(move.processed,'exclude')
                excludem(find(move.processed.exclude.SFPs) == 0) = 0;
                outofFOVm.exclude_badSFP = move.processed.exclude.SFPs;
            end
        end
%         allPrepped{si,1} = allPrepped{si,1}(excluder,:,:);
        allPrepped = allPrepped(excludem,:,:);
        outofFOVr.exclude_all{si} = excluder;
        outofFOVm.exclude_all{si} = excludem;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        itermap = repmat({[]},[1 iters]);
        iterstruct = repmat({[]},[1 iters]);
        
        if combs(si,1) > 9
            preppedr = permute(sessions{combs(si,1)}.SFPs,[3 1 2]);
            preppedr = preppedr(outofFOVr.exclude_all{si},:,:);
            outP = [path '\Segments\ms' num2str(combs(si,1))];% num2str(find(si == combs(i,:)))];
            checkP(outP)
            save(outP,'preppedr');
        else
            preppedr = permute(sessions{combs(si,1)}.SFPs,[3 1 2]);
            preppedr = preppedr(outofFOVr.exclude_all{si},:,:);
            outP = [path '\Segments\ms0' num2str(combs(si,1))];% num2str(find(si == combs(i,:)))];
            checkP(outP)
            save(outP,'preppedr');
        end
        if combs(si,2) > 9            
            preppedm = allPrepped;
            outP = [path '\Segments\ms' num2str(combs(si,2))];% num2str(find(si == combs(i,:)))];
            checkP(outP)
            save(outP,'preppedm');
        else
            preppedm = allPrepped;
            outP = [path '\Segments\ms0' num2str(combs(si,2))];% num2str(find(si == combs(i,:)))];
            checkP(outP)
            save(outP,'preppedm');
        end
        
        mkdir(path,'tempfigs')
        parfor iteration = 1:iters
            try
            if nonRigid
                [itter_map, regStruct] = registerCellsNoGUI2020_Auto([path '\Segments' ]);
            else
                [itter_map, regStruct] = registerCells2020_Auto([path '\Segments' ]);
            end
            catch            
                itter_map = [];
                regStruct = [];
            end
            itermap{iteration} = itter_map;
            iterstruct{iteration} = regStruct;
            close all
            close all hidden
            drawnow
            copyfile([path,'\Segments\Plots\CellRegistration\Figures'],[path,'\tempfigs\',num2str(iteration)]);
        end
        try
            rmdir([path, '\Segments' ],'s');
        end        
        iterstruct = iterstruct(~cellfun(@isempty,itermap));
        itermap = itermap(~cellfun(@isempty,itermap));
        itterind = find(~cellfun(@isempty,itermap));
        
        [a, b] = nanmin(cellfun(@length,itermap));
        %             [a b] = nanmax(scores);
        if ~isempty(b)
            map = itermap{b};
            score = iterstruct{b}.cell_scores;
            corval = iterstruct{b}.maximal_cross_correlation;
            if ~nonRigid
                corval_foot = iterstruct{b}.maximal_cross_correlation_foot;
            end
            for s = 1 : length(iterstruct{b}.p_same_registered_pairs)
                probtemp(s,:) = nanmean(iterstruct{b}.p_same_registered_pairs{s});
            end
            copyfile([path,'\tempfigs\',num2str((itterind(b)))],[path,'\Results\',num2str(combs(si,1)),'_',num2str(combs(si,2))]);
            %                 corval = iterstruct{b}.maximal_cross_correlation;
        end
        
        alignmentMap{si} = map;
        prob{si} = probtemp;
        scoreMap{combs(si,1),combs(si,2)} = score;
        cormat(si) = corval;
        if ~nonRigid
            cormat_foot = corval_foot;
        end
        probtemp = [];
        try
            rmdir([path, '\tempfigs' ],'s');
        end
        
        best = [alignmentMap(si)];% oldAlignmentMap(combs(i,1),combs(i,2))];
        best = best(~cellfun(@isempty,best));
        bprob = [prob(si)];% oldAlignmentMap(combs(i,1),combs(i,2))];
        bprob = bprob(~cellfun(@isempty,bprob));
        bestscore = [scoreMap(combs(si,1),combs(si,2))];% oldAlignmentMap(combs(i,1),combs(i,2))];
        bcorr = [cormat(si)];% oldAlignmentMap(combs(i,1),combs(i,2))];
        if ~nonRigid
            bcorr_foot = [cormat_foot];
        end
        bestscore = bestscore(~cellfun(@isempty,bestscore));
        if isempty(best)
            continue
        end
        [a b] = nanmin(cellfun(@length,best));
        nm = best;
        nmp = bprob;
        nms = bestscore;
        nmc_cent = bcorr;
        if ~nonRigid
            nmc_foot = bcorr_foot;
        end
        alignmentMap(si) = nm;
        probMap(si) = nmp;
        scoreMap(combs(si,1),combs(si,2)) = nmp;
        CorrMap_cent(combs(si,1),combs(si,2)) = nmc_cent;
        if ~nonRigid
            CorrMap_foot(combs(si,1),combs(si,2)) = nmc_foot;
        end
    end
    %{
    for i = 1:length(combs(:,1))        
        itermap = repmat({[]},[1 iters]);
        iterstruct = repmat({[]},[1 iters]);
        
        preppedr = permute(sessions{combs(i,1)}.SFPs,[3 1 2]);
        preppedr = preppedr(outofFOVr.exclude_all{i},:,:);
        outP = [path '\Segments\ms' num2str(combs(i,1))];% num2str(find(si == combs(i,:)))];
        checkP(outP)
        save(outP,'preppedr');
        
        preppedm = allPrepped{i};
        outP = [path '\Segments\ms' num2str(combs(i,2))];% num2str(find(si == combs(i,:)))];
        checkP(outP)
        save(outP,'preppedm');
        
        mkdir(path,'tempfigs')
        for iteration = 1:iters
            %                 try
            if nonRigid
                [map, regStruct] = registerCellsNoGUI2020([path '\Segments' ]);
            else
                [map, regStruct] = registerCells2020([path '\Segments' ]);
            end
            %                 catch
            %                     map = [];
            %                     regStruct = [];
            %                 end
            itermap{iteration} = map;
            iterstruct{iteration} = regStruct;
            close all
            close all hidden
            drawnow
            copyfile([path,'\Segments\Plots\CellRegistration\Figures'],[path,'\tempfigs\',num2str(iteration)]);
        end
        try
            rmdir([path, '\Segments' ],'s');
        end
        iterstruct = iterstruct(~cellfun(@isempty,itermap));
        itermap = itermap(~cellfun(@isempty,itermap));
        
        [a, b] = nanmin(cellfun(@length,itermap));
        %             [a b] = nanmax(scores);
        if ~isempty(b)
            map = itermap{b};
            score = iterstruct{b}.cell_scores;
            corval = iterstruct{b}.maximal_cross_correlation;
            if ~nonRigid
                corval_foot = iterstruct{b}.maximal_cross_correlation_foot;
            end
            for s = 1 : length(iterstruct{b}.p_same_registered_pairs)
                probtemp(s,:) = nanmean(iterstruct{b}.p_same_registered_pairs{s});
            end
            copyfile([path,'\tempfigs\',num2str(b)],[path,'\Results\',num2str(combs(i,1)),'_',num2str(combs(i,2))]);
            %                 corval = iterstruct{b}.maximal_cross_correlation;
        end
        alignmentMap{i} = map;
        prob{i} = probtemp;
        scoreMap{i} = score;
        cormat(i) = corval;
        if ~nonRigid
            cormat_foot(i) = corval_foot;
        end
        probtemp = [];
        try
            rmdir([path, '\tempfigs' ],'s');
        end
    end
    if nFold == 2
        nm = repmat({[]},length(sessions));
        nmp = repmat({[]},length(sessions));
        nms = repmat({[]},length(sessions));
        nmc_cent = repmat([],length(sessions));
        nmc_foot = repmat([],length(sessions));
        for i = 1:length(combs(:,1))
            best = [alignmentMap(i)];% oldAlignmentMap(combs(i,1),combs(i,2))];
            best = best(~cellfun(@isempty,best));
            bprob = [prob(i)];% oldAlignmentMap(combs(i,1),combs(i,2))];
            bprob = bprob(~cellfun(@isempty,bprob));
            bestscore = [scoreMap(i)];% oldAlignmentMap(combs(i,1),combs(i,2))];
            bcorr = [cormat(i)];% oldAlignmentMap(combs(i,1),combs(i,2))];
            if ~nonRigid
                bcorr_foot = [cormat_foot(i)];
            end
            bestscore = bestscore(~cellfun(@isempty,bestscore));
            if isempty(best)
                continue
            end
            [a b] = nanmin(cellfun(@length,best));
            nm{combs(i,1),combs(i,2)} = best{b};
            nmp{combs(i,1),combs(i,2)} = bprob{b};
            nms{combs(i,1),combs(i,2)} = bestscore{b};
            nmc_cent(combs(i,1),combs(i,2)) = bcorr(b);
            if ~nonRigid
                nmc_foot(combs(i,1),combs(i,2)) = bcorr_foot(b);
            end
        end
        alignmentMap = nm;
        probMap = nmp;
        scoreMap = nms;
        CorrMap_cent = nmc_cent;
        if ~nonRigid
            CorrMap_foot = nmc_foot;
        end
    else
        nm = repmat({[]},[length(combs(:,1)) 1]);
        nmp = repmat({[]},[length(combs(:,1)) 1]);
        nms = repmat({[]},[length(combs(:,1)) 1]);
        nmc_cent = repmat([],[length(combs(:,1)) 1]);
        if ~nonRigid
            nmc_foot = repmat([],[length(combs(:,1)) 1]);
        end
        for i = 1:length(combs(:,1))
            best = [alignmentMap(i)];% oldAlignmentMap(i)];
            best = best(~cellfun(@isempty,best));
            bprob = [prob(i)];% oldAlignmentMap(i)];
            bprob = bprob(~cellfun(@isempty,bprob));
            bestscore = [scoreMap(i)];% oldAlignmentMap(i)];
            bestscore = bestscore(~cellfun(@isempty,bestscore));
            bcorr = [cormat(i)];% oldAlignmentMap(i)];
            if ~nonRigid
                bcorr_foot = [cormat_foot(i)];
            end
            if isempty(best)
                continue
            end
            [a b] = nanmin(cellfun(@length,best));
            nm{i} = best{b};
            nmp{i} = bprob{b};
            nms{i} = bestscore{b};
            nmc_cent(i) = bcorr(b);
            if ~nonRigid
                nmc_foot(i) = bcorr_foot(b);
            end
        end
        alignmentMap = nm;
        probMap = nmp;
        scoreMap = nms;
        CorrMap_cent = nmc_cent;
        if ~nonRigid
            CorrMap_foot = nmc_foot;
        end
    end
    %}
    
            
            
    alignmentMap = updateAlignment(alignmentMap,outofFOVm,combs);
    Singlemap = ReorganizeAlignmentMap(alignmentMap);
    [Singlemap, avg_psame] = ElimConflict2020(Singlemap,alignmentMap,scoreMap);
    Singlemap = FindMissingCells2020_Auto(Singlemap,outofFOVm,alignmentMap,combs);
    
    alignment.alignmentMap = alignmentMap;
    alignment.Singlemap = Singlemap;
    alignment.combs = combs;
    alignment.sessions = sessions;
    alignment.nFold = nFold;
    alignment.probMap = probMap;
    alignment.scoreMap = scoreMap;
    if ~nonRigid
        alignment.CorrMap_footprints = CorrMap_foot;
        alignment.CorrMap_centroid = CorrMap_cent;
    else
        alignment.CorrMap = CorrMap_cent;
    end
    %         save('temporaryAlignmentMap','alignmentMap')
    save('alignment','-struct','alignment','-v7.3');
% end
end