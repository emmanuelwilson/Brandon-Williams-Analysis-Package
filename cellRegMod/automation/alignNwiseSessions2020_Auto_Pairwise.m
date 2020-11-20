%% Modified cellReg script, combination of Alex's and Emmmanuel's

function alignNwiseSessions2020_Auto_Pairwise(SFPsShifted,sessions,combs,path,nonRigid)
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
        load(['ms' num2str(combs(si,1)) '.mat'])
        ref = ms;
        preppedr = ref.SFPs .* bsxfun(@gt,ref.SFPs,0.5.*nanmax(nanmax(ref.SFPs,[],1),[],2));
        move = SFPsShifted{si};%,'calcium','processed');        
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
            preppedr = permute(ref.SFPs,[3 1 2]);
            preppedr = preppedr(outofFOVr.exclude_all{si},:,:);
            outP = [path '\Segments\ms' num2str(combs(si,1))];% num2str(find(si == combs(i,:)))];
            checkP(outP)
            save(outP,'preppedr');
        else
            preppedr = permute(ref.SFPs,[3 1 2]);
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
%             try
            if nonRigid
                [itter_map, regStruct] = registerCellsNoGUI2020_Auto([path '\Segments' ]);
            else
                [itter_map, regStruct] = registerCells2020_Auto([path '\Segments' ]);
            end
%             catch            
%                 itter_map = [];
%                 regStruct = [];
%             end
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
            while isempty(itermap{b})
                itermap(b) = [];
                iterstruct(b) = [];
                itterind(b) = [];
                [a, b] = nanmin(cellfun(@length,itermap));
            end
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
        else 
            nothing = 0;
        end
        
        alignmentMap{si} = map;
        prob{si} = probtemp;
        scoreMap{si,1} = score;
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
        bestscore = [scoreMap(si,1)];% oldAlignmentMap(combs(i,1),combs(i,2))];
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
        try
            probMap(si) = nmp;
        end
        scoreMap(si,1) = nmp;
        CorrMap_cent(si,1) = nmc_cent;
        if ~nonRigid
            CorrMap_foot(si,1) = nmc_foot;
        end
    end
    
    alignmentMap = updateAlignment(alignmentMap,outofFOVm,combs);
    
    alignmentMap = alignmentMap;    
    combs = combs;   
    nFold = nFold;
    probMap = probMap;
    scoreMap = scoreMap;
    if ~nonRigid
        alignment.CorrMap_footprints = CorrMap_foot;
        alignment.CorrMap_centroid = CorrMap_cent;
    else
        alignment.CorrMap = CorrMap_cent;
    end
    %         save('temporaryAlignmentMap','alignmentMap')
    save('alignment','alignment','combs','nFold','probMap','scoreMap','-v7.3');
% end
end