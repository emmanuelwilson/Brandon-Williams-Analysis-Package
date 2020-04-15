function mapAnalysis(paths,doPlot)

    if nargin < 2 || isempty(doPlot)
        doPlot = false;
    end

    clc
    fprintf('\n')
    %%% Reliability constaint
    
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
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
    
    labels = [{'Saline'} {'CNO'} {'No Injection'}];
    velThresh = -2;
    pThresh = 0.05;
    afw = repmat({[]},[length(labels) length(upiece)]);
    amfr = repmat({[]},[length(labels) length(upiece)]);
    amfr = repmat({[]},[length(labels) length(upiece)]);
    icorr = repmat({[]},[length(labels) length(upiece)]);
    cellCounts = repmat({[]},[length(labels) length(upiece)]);
    crossAllComp = repmat({[]},[length(labels) length(upiece)]);
    for mi = 1:length(upiece)
        fprintf(['\n\n\tMouse:  ' num2str(upiece{mi}) '\n']) 
        isM = find(ismember(piece,upiece(mi)));
        for si = 1:length(isM);
            fprintf(['\n\t\tSession:  ' paths{isM(si)}])
            s = load(paths{isM(si)});
%             s.processed.p = imfilter(s.processed.p,fspecial('gauss',[1 15],2),'same','replicate');
            %id group
            group = find(ismember(labels,paths{isM(si)}(find(ismember(paths{isM(si)},'_'),1,'last')+1:end-4)));
            
            if isempty(group)
                group = 3;
            end
            
            %compute velocity
            vel = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
            vel = imfilter(vel,fspecial('gauss',[1 30],10),'same','replicate');

            cellType = nan(1,4);
            
            %choose cells
            if isfield(s.processed,'exclude')
                inds = (s.processed.splithalf.p<=pThresh & s.processed.exclude.SFPs);
                inc = s.processed.exclude.SFPs;
                cellType = [nansum(inc & s.processed.splithalf.p>pThresh) ...
                    nansum(inc & s.processed.splithalf.p<=pThresh)];
            else
                
            end
            
            cellCounts{group,mi,1} = [cellCounts{group,mi,1}; {cellType(1)} {cellType(2)}];
            str = sprintf(['\n\t\t\tCell Breakdown:  %0.2f%%  %0.2f%%  %0.2f%%  %0.2f%%  '],100.*cellType./nansum(cellType));
            fprintf('%s',str);
            
             gT = s.processed.trace(inds,:);
            
            partitions = [1; find(diff(s.processed.validTraceFrames(:,2))>300); length(s.processed.validTraceFrames(:,2))+1];
            
            for pi = 1:length(partitions)-1
                mask = false(1,length(s.processed.p(1,:)));
                mask(partitions(pi):partitions(pi+1)-1) = true;
                allMasks{pi} = [mask & vel>velThresh];
            end
            
            [map samp allComp] = getMatchedMapsNMasks(s.processed.p,gT,allMasks);

            crossAllComp{group,mi} = cat(1,crossAllComp{group,mi},...
                [help_getMaskedVals(allComp(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(allComp(1:2,3:4,:),[false true; true false])]);
            
            if doPlot
                map = getMatchedMapsNMasks(s.processed.p,gT,allMasks);
                help_plotMaps(map,paths(isM(si)));
            end
            
            [blah1 ivals mfr pfr] = getMatchedSamplingValues(s.processed.p,gT,allMasks);
            tmp = [help_getMaskedVals(mfr(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(mfr(1:2,3:4,:),[false true; true false])];
            amfr(group,mi) = {[amfr{group,mi}; tmp]};
            icorr(group,mi) = {[icorr{group,mi}; help_getMaskedVals(ivals(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(ivals(1:2,3:4,:),[false true; true false])]};
            
        end
    end
    close all
    
    eliminate = all(cellfun(@isempty,crossAllComp),2);
    labels(eliminate) = [];
    amfr(eliminate,:) = [];
    icorr(eliminate,:) = [];
    afw(eliminate,:) = [];
    cellCounts(eliminate,:) = [];
    crossAllComp(eliminate,:) = [];
    
    
    slashInds = find(ismember(paths{1},'/'));
    root = ['Plots/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
    
    h1 = mkPie(cellCounts);
%     lgd = legend([h1 h2],[{'Nonspatial'} {'Spatial'} {'Room'} {'Hallway'} {'Both'}],...
%         'box','off','location','northoutside','orientation','horizontal');
    saveFig(gcf,[root '/CellClassification'],[{'pdf'} {'tiff'}]);

    figure
    set(gcf,'position',[50 450 175.*length(upiece) 200])
    toPlot = repmat({[]},[2 length(labels) length(upiece)]);
    lim = cat(1,crossAllComp{:});
    lim = ceil(nanmax(lim(:)).*10)./10;
    for ui = 1:length(upiece)
        for group = 1:length(labels)
            toPlot{1,group,ui} = crossAllComp{group,ui}(:,1);
            toPlot{2,group,ui} = crossAllComp{group,ui}(:,2);
        end
        subplot(1,length(upiece),ui)
        mkGraph(toPlot(:,:,ui))
        ylabel('Correlation (r)')
        tmp = upiece{ui};
        slashInds = find(ismember(tmp,'/'));
        set(gca,'xticklabel',{upiece{ui}(slashInds(end)+1:end)})
        set(gca,'ylim',[0 lim])
    end
    saveFig(gcf,[root '/PopVecRemapping_AnimalWise'],[{'pdf'} {'tiff'}]);
    
    figure
    set(gcf,'position',[50 450 175 200])
    toPlot = repmat({[]},[2 length(labels)]);
    diffPlot = repmat({[]},[1 length(labels)]);
    for group = 1:length(labels)
        tmp = cat(1,crossAllComp{group,:});
        toPlot{1,group} = tmp(:,1);
        toPlot{2,group} = tmp(:,2);
        diffPlot{group} = tmp(:,1)-tmp(:,2);
    end
    mkGraph(toPlot)
    ylabel('Correlation (r)')
    set(gca,'ylim',[0 lim])
    saveFig(gcf,[root '/PopVecRemapping_Combined'],[{'pdf'} {'tiff'}]);
    
    figure
    set(gcf,'position',[50 450 175 200])
    mkGraph(diffPlot)
    ylabel('Correlation (r)')
%     set(gca,'ylim',[0 lim])
    saveFig(gcf,[root '/PopVecRemapping_Combined_Difference'],[{'pdf'} {'tiff'}]);
        
    figure
    set(gcf,'position',[500 450 300.*length(upiece) length(labels).*300])
    for ui = 1:length(upiece)
        for i = 1:length(labels)
            subplot(length(labels),length(upiece),(i-1).*length(upiece)+ui)
            h = cumHist(amfr{i,ui},[0:0.025:1]);
            xlabel('Mean Rate Change (%)')
            ylabel('Count')
            legend([h{1}(1) h{2}(1)],[{'Within'} ...
                {'Across'} {'Shuffled'}],'location','southeast','fontname','arial',...
                'fontsize',9,'fontweight','bold','color','none','box','off')
            set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
            tmp = upiece{ui};
            slashInds = find(ismember(tmp,'/'));
            title({upiece{ui}(slashInds(end)+1:end)})
            axis square
        end
    end
    drawnow
    saveFig(gcf,[root '/RateRemapping_AnimalWise'],[{'tiff'} {'pdf'}]);
    
    figure
    set(gcf,'position',[500 450 300 length(labels).*300])
    for i = 1:length(labels)
        subplot(length(labels),1,i)
        h = cumHist(cat(1,amfr{i,:}),[0:0.025:1]);
        xlabel('Mean Rate Change (%)')
        ylabel('Count')
        legend([h{1}(1) h{2}(1)],[{'Within'} ...
            {'Across'} {'Shuffled'}],'location','southeast','fontname','arial',...
            'fontsize',9,'fontweight','bold','color','none','box','off')
        set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
        axis square
    end
    drawnow
    saveFig(gcf,[root '/RateRemapping_Combined'],[{'tiff'} {'pdf'}]);

    
    outP = ['Stats/' paths{1}(slashInds(1):slashInds(2)-1) '.txt'];
    checkP(outP);
    fid = fopen(outP,'w');
    fprintf(fid,'\t\t\tPop Vec Remapping t-test\n');
    allGroupComps = [];
    for group = 1:length(labels)
        tmp = cat(1,crossAllComp{group,:});
        fprintf(fid,['\nWithin Group:  ' labels{group}]);
        [h pval ci tstat] = ttest(tmp(:,1),tmp(:,2));
        fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.9f ',tstat.df,tstat.tstat,pval);
        allGroupComps = [allGroupComps {diff(fliplr(tmp),[],2)}];
    end
    for i = 1:length(allGroupComps)
        for j = i+1:length(allGroupComps)
            fprintf(fid,['\nAcross Group:  ' labels{i} ' vs. ' labels{j}]);
            [h pval ci tstat] = ttest2(allGroupComps{i},allGroupComps{j});
            fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.9f ',tstat.df,tstat.tstat,pval);
        end
    end
    
    tmp = cat(1,amfr{i,:});
    fprintf(fid,'\t\t\tCell-wise rate change t-test\n');
    allGroupComps = [];
    for group = 1:length(labels)
        tmp = cat(1,amfr{i,:})
        fprintf(fid,['\nWithin Group:  ' labels{group}]);
        [h pval ci tstat] = ttest(tmp(:,1),tmp(:,2));
        fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.9f ',tstat.df,tstat.tstat,pval);
        allGroupComps = [allGroupComps {diff(fliplr(tmp),[],2)}];
    end
    
    fclose(fid);
    
    out.labels = labels;
    out.amfr = amfr;
    out.icorr = icorr;
    out.afw = afw;
    out.cellCounts = cellCounts;
    out.crossAllComp = crossAllComp;
    
    dataP = ['AggregatedData' paths{1}(slashInds(1):slashInds(2)-1)];
    checkP(dataP);
    save(dataP,'-struct','out','-v7.3');
end

function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end

function help_plotMaps(totalMaps,p)

    doK = [8 4];

    for part = 0:floor(length(totalMaps(1,1,:,1))/prod(doK))

        figure(1)
        set(gcf,'position',[50 50 900 1350])
        for k = 1:prod(doK)
            if part.*prod(doK)+k > length(totalMaps(1,1,:,1))
                break
            end
%             tmp = [totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,2)];

            tmp = [totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,3); ...
                nan(1,length(totalMaps(1,:,1,1)).*2+1); ...
                totalMaps(:,:,part.*prod(doK)+k,2) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,4)];

            subplot(doK(1),doK(2),k)

            imagesc(tmp)
            colormap jet
            caxis([0 nanmax(tmp(:))])
            alpha(double(~isnan(tmp)))
            axis equal
            axis off    
        end

        slashInds = find(ismember(p{1},'/'));
        outP = ['Plots/DifferentiatedMaps/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
        saveFig(gcf,outP,'tiff')
        saveFig(gcf,outP,'pdf')
        close all
        drawnow
    end
end
























