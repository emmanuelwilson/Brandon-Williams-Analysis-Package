function pairwiseMapAnalysis(paths,doPlot)
    
    close all
    drawnow
    
    params = ebcMapParams();
    
    pThresh = 0.05; % 0.025
    
    if nargin < 2
        doPlot = false;
    end

    warning off all
    
    clc
    fprintf('\nComputing pairwise map analysis ')
    if doPlot
        fprintf('and plotting maps')
    end
    fprintf('\n')
    %%% Reliability constaint
    
    %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    labels = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
        blah = find(ismember(piece{end},'_'),2,'last');
    end
    upiece = unique(piece);
    
    mouseName = upiece;
    for mi = 1:length(upiece)
        mouseName{mi} = upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end);
    end
    
    envLabel = [{'A'} {'A2'} {'B'} {'A3'}];
    doComp = [{'A'} {'A2'}; {'A2'} {'B'}; {'B'} {'A3'}; ...
        {'A'} {'B'}; {'A2'} {'A3'}; {'A'} {'A3'}]; %%% Complete curve

%     doComp = [{'A'} {'A2'}; {'A'} {'A3'}; {'A2'} {'A3'}; ...
%         {'A'} {'B'}; {'A2'} {'B'}; {'A3'} {'B'}]; %%% Complete curve
    
    aapfd = [];
    aad = [];
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        doAl = help_getAlignmentID(s.alignment,2,paths(isM));
        if isnan(doAl)
            fprintf(['\n\t\t***** No pairwise alignment for these sessions *****\n'])
            continue
        end
        tic
        alignMap = s.alignment(doAl(1)).alignmentMap;
        envs = [];
        am = repmat({[]},[1 length(sessions)]);
        apfd = repmat({[]},[1 length(sessions)]);
        ad = repmat({[]},[1 length(sessions)]);
        awm = repmat({[]},[1 length(sessions)]);
        isPC = repmat({[]},[1 length(sessions)]);
        fprintf(['\t\tLoading maps... '])
        sequenceNum = str2num(upiece{mi}(find(ismember(upiece{mi},'_'),1,'last')+1:end));
        for si = 1:length(sessions)
            s = load(sessions{si},'processed','properties','exclude');
            if isfield(s.processed,'exclude')
                exclude = s.processed.exclude.SFPs;
            else
                exclude = true(length(s.processed.splithalf.wholeEBC),1);
            end
            am{si} = s.processed.ebc.whole;
            awm{si} = s.processed.ebc.walls;
            
%             doInclude = s.processed.splithalf.ebc.mrl.p < pThresh & ...
%                 s.processed.splithalf.ebc.shc.p < pThresh;
            
            doInclude = s.processed.splithalf.ebc.mrl.p < 0.05 & ...
                s.processed.ebc.shPFD < 45 & s.processed.ebc.shD < 3; % & s.processed.ebc.shD < 3
                
            isPC{si} = doInclude & exclude;
            
            apfd{si} = s.processed.ebc.pfd(doInclude & exclude);
            ad{si} = s.processed.ebc.mrld(doInclude & exclude);
            
            envs = [envs; {s.properties.session}];
        end  
        
        durat = toc;
        fprintf([num2str(durat) ' s\n']);
        
        for si = 1:length(sessions)
            for sj = 1:length(sessions)
                if si == sj
                    alignMap{si,sj} = [[1:length(am{si}(1,1,:))]' [1:length(am{si}(1,1,:))]'];
                end                
                if isempty(alignMap{si,sj})
                    alignMap{si,sj} = alignMap{sj,si}(:,[2 1]);
                end
                if isempty(alignMap{si,sj})
                    continue
                end
                alignMap{si,sj} = alignMap{si,sj}(all(alignMap{si,sj}~=0,2),:);
            end
        end
        
        tic
        iter = 0;
        strLength = 0;
        fprintf(['\t\tComputing map comparisons...\n\t\t\t'])
        crossSim = nan(length(sessions));
        kCrossSim = repmat({[]},length(sessions));
        kCrossSim_shuffle = repmat({[]},length(sessions));
        kAngDiff = repmat({[]},length(sessions));
        kAngDiff_shuffle = repmat({[]},length(sessions));
        kMaxCorr = repmat({[]},length(sessions));
        kMaxCorr_shuffle = repmat({[]},length(sessions));
        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                
                iter = iter+1;
                fprintf(repmat('\b',[1 strLength]));
                str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2))]);
                
                fprintf(str);
                strLength = length(str);
                
                if isempty(alignMap{si,sj})
                    continue
                end
                
                isGood = isPC{si}(alignMap{si,sj}(:,1)) | isPC{sj}(alignMap{si,sj}(:,2));

                tmp1 = am{si}(:,:,alignMap{si,sj}(isGood,1));
                tmp2 = am{sj}(:,:,alignMap{si,sj}(isGood,2));
               
                rt1 = reshape(tmp1,[numel(tmp1(:,:,1)) length(tmp1(1,1,:))]);
                rt2 = reshape(tmp2,[numel(tmp2(:,:,1)) length(tmp2(1,1,:))]);
                
                isGP = ~[any(isnan(rt1),2)|any(isnan(rt2),2)];
                
                xc = corr(rt1(isGP,:),rt2(isGP,:));
                kCrossSim{si,sj} = xc(logical(eye(size(xc))));
                
                [kAngDiff{si,sj} kMaxCorr{si,sj}] = getAngDiff(tmp1,tmp2);
                
                
                %%% Shuffled sims
                nsim = 100;
                null = nan(length(rt2(1,:)),nsim);
                null_ad = nan(length(rt2(1,:)),nsim);
                null_maxCorr = nan(length(rt2(1,:)),nsim);
                for i = 1:nsim
                    rt2 = rt2(:,randperm(length(rt2(1,:))));
                    xc = corr(rt1(isGP,:),rt2(isGP,:));
                    null(:,i) = xc(logical(eye(size(xc))));
                    
                    [null_ad(:,i) null_maxCorr(:,i)] =  ...
                        getAngDiff(tmp1,tmp2(:,:,randperm(length(tmp2(1,1,:)))));
                end
                kCrossSim_shuffle{si,sj} = null(:);
                kAngDiff_shuffle{si,sj} = null_ad(:);
                kMaxCorr_shuffle{si,sj} = null_maxCorr(:);
            end
        end
        
        durat = toc;
        fprintf(['\t' num2str(durat) ' s\n']);
        
        root = ['Plots/Summary/' ...
            upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)];
        
        figure
        transitionPlot(kMaxCorr,envs,doComp);
        h = transitionPlot(kMaxCorr_shuffle,envs,doComp,false);
        set(h,'color',[0.5 0.5 0.5],'linestyle','--')
        ylabel('Corrected EBC Map Corr (r)');
        outP = [root '/MaxMapCorr'];
        saveFig(gcf,outP,[{'tiff'} {'pdf'}])
        
%         figure
%         transitionPlot(kCrossSim,envs,doComp);
%         h = transitionPlot(kCrossSim_shuffle,envs,doComp,false);
%         set(h,'color',[0.5 0.5 0.5],'linestyle','--')
        
        figure
        transitionPlot(kAngDiff,envs,doComp);
        h = transitionPlot(kAngDiff_shuffle,envs,doComp,false);
        set(h,'color',[0.5 0.5 0.5],'linestyle','--')
        ylabel('Orientation Difference (deg)');
        outP = [root '/Orientation'];
        saveFig(gcf,outP,[{'tiff'} {'pdf'}])  
        
        [a b] = ismember(envs,envLabel);
        aapfd = [aapfd; apfd(b)];
        aad = [aad; ad(b)];
        close all
    end
    
    figure
    h = mkGraph(aapfd,envLabel);
    legend(h,mouseName','location','southoutside','orientation','horizontal')
    set(h,'linestyle','none')
    set(gca,'ylim',[-30 390],'ytick',[0:30:360])
    hold on
    for i = 1:5
        plot(get(gca,'xlim'),(i-1).*90.*ones(1,2),'linestyle','--', ...
            'color',[0.5 0.5 0.5])
    end
    root = ['Plots/Summary/PFD_Distribution'];
    
    
    figure
    h = mkGraph(aad,envLabel);
    legend(h,mouseName','location','southoutside','orientation','horizontal')
    set(h,'linestyle','none')
    set(gca,'ylim',[-30 390],'ytick',[0:30:360])
    hold on
    for i = 1:5
        plot(get(gca,'xlim'),(i-1).*90.*ones(1,2),'linestyle','--', ...
            'color',[0.5 0.5 0.5])
    end
    root = ['Plots/Summary/PFD_Distribution'];
end


















