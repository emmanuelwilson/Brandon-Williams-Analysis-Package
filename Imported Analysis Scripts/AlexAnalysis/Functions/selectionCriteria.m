function selectionCriteria(paths,varargin)
    clc
    close all
    drawnow
    
    if isempty(varargin)
        varargin = [{'ebc_shc'}];
    end
    varargin = cellfun(@lower,varargin,'uniformoutput',false);
    
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
    
    fprintf(['Computing split-half reliability:\n']);
    tmp = [repmat({'\n\t'},[1 length(varargin)]); varargin];
    fprintf(cat(2,tmp{:},'\n\n'))
    velThresh = -2;
    nsims = 1000; %500
    minShift = 900;
    for p = paths'
        s = load(p{1});
        didChange = false;

        fprintf(['\t' num2str(p{1}) '\n'])    
        if isfield(s.processed,'exclude')
            include = s.processed.exclude.SFPs;
        else
            include = true(length(s.processed.trace(:,1)),1);
        end
        
        if ismember({'ebc_shc'},varargin)

            didChange = true;
            
            fprintf('\t\tWhole Map (Unmatched; EBC), computing null...  ')
            
            tic
            [pval actual null] = getshEBC(s.processed.p,s.processed.hd,s.processed.trace);
            
            s.processed.splithalf.ebc.shc.val = actual;
            s.processed.splithalf.ebc.shc.p = pval;
            s.processed.splithalf.ebc.shc.null = null;
            
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            
            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.ebc.shc.p(include),[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/WholeEBCMap/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.ebc.shc.p(include) <= 0.05) ...
                length(s.processed.splithalf.ebc.shc.p(include)) nanmean(s.processed.splithalf.ebc.shc.p(include) <= 0.05).*100]);
        end
        
        if ismember({'ebc_mrl'},varargin)

            didChange = true;
            
            allMasks = [{true(1,length(s.processed.p(1,:)))}];
            fprintf('\t\tWhole Map (Unmatched; EBC), computing null...  ')
            tic
            [pval actual null] = getUnmatchedEBCNMasks( ...
                s.processed.p,s.processed.hd,s.processed.trace,allMasks);
            
            s.processed.splithalf.ebc.mrl.val = actual;
            s.processed.splithalf.ebc.mrl.p = pval;
            s.processed.splithalf.ebc.mrl.null = null;
            
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            
            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.ebc.mrl.p(include),[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/WholeEBCMap/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.ebc.mrl.p(include) <= 0.05) ...
                length(s.processed.splithalf.ebc.mrl.p(include)) nanmean(s.processed.splithalf.ebc.mrl.p(include) <= 0.05).*100]);
        end
        
        if ismember({'points_mrl'},varargin)

            didChange = true;
            
            fprintf('\t\tWhole Map (Unmatched; EBC), computing null...  ')
            tic
            [pval actual null] = getShuffledPointMaps( ...
                s.processed.p,s.processed.hd,s.processed.trace,s.processed.locs);
            
            s.processed.splithalf.points.mrl.val = actual;
            s.processed.splithalf.points.mrl.p = pval;
            s.processed.splithalf.points.mrl.null = null;
            
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            
            figure(1)
            set(gcf,'position',[50 50 250 250])
            toPlot = repmat({[]},[1 length(s.processed.splithalf.points.mrl.p(1,:))]);
            for i = 1:length(s.processed.splithalf.points.mrl.p(1,:))
                toPlot{i} = s.processed.splithalf.points.mrl.p(include,i);
            end
            cumHist(toPlot,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/WholePointMaps/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.points.mrl.p(include,1) <= 0.05) ...
                length(s.processed.splithalf.points.mrl.p(include,1)) nanmean(s.processed.splithalf.points.mrl.p(include,1) <= 0.05).*100]);
        end
        
        
        
        if ismember({'wholemap_si'},varargin)

            didChange = true;
            
            allMasks = [{true(1,length(s.processed.p(1,:)))}];
            fprintf('\t\tWhole Map (Unmatched; Spatial Information), computing null...  ')
            tic
            [pval actual null] = getUnmatchedSICircNullNMasks(s.processed.p,s.processed.trace,allMasks);
            
            s.processed.splithalf.wholemap_unmatched_si.val = actual;
            s.processed.splithalf.wholemap_unmatched_si.p = pval;
            s.processed.splithalf.wholemap_unmatched_si.null = null;
            
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            
            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.wholemap_unmatched_si.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/WholeMap_SI/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.wholemap_unmatched_si.p(include) <= 0.05) ...
                length(s.processed.splithalf.wholemap_unmatched_si.p(include)) nanmean(s.processed.splithalf.wholemap_unmatched_si.p(include) <= 0.05).*100]);
        end
        
        if ismember({'wholemap_unmatched'},varargin)
            %%% In-room split half reliability, ignore doorway

%             if isfield(s.processed,'splithalf') && isfield(s.processed.splithalf,'wholemap')
%                 continue
%             end
            didChange = true;
            
            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
            allMasks = [{half} {~half}];
            fprintf('\t\tWhole Map (Unmatched), computing null...  ')
            tic
            [pval actual null] = getUnmatchedCircNullNMasks(s.processed.p,s.processed.trace,allMasks);
            
            s.processed.splithalf.wholemap_unmatched.val = permute(actual(1,2,:),[3 1 2]);
            s.processed.splithalf.wholemap_unmatched.p = permute(pval(1,2,:),[3 1 2]);
            s.processed.splithalf.wholemap_unmatched.null = permute(null(1,2,:),[3 1 2]);
            
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            
            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.wholemap_unmatched.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/WholeMap_Unmatched/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.wholemap_unmatched.p(include) <= 0.05) ...
                length(s.processed.splithalf.wholemap_unmatched.p(include)) nanmean(s.processed.splithalf.wholemap_unmatched.p(include) <= 0.05).*100]);
        end

        if ismember({'wholemap'},varargin)
            %%% In-room split half reliability, ignore doorway

%             if isfield(s.processed,'splithalf') && isfield(s.processed.splithalf,'wholemap')
%                 continue
%             end
            didChange = true;
            
            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
            allMasks = [{half} {~half}];
            [a b c ival] = getMatchedMapsNMasks(s.processed.p,s.processed.trace,allMasks);
            s.processed.splithalf.wholemap.val = ival;

            null = nan(length(s.processed.trace(:,1)),nsims);
            fprintf('\t\tWhole Map, computing null...  ')
            tic
            P = s.processed.p;
            T = s.processed.trace;
            parfor sim = 1:nsims
                gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);

                [map samp allComp ival] = getMatchedMapsNMasks(P,gT,allMasks);
                null(:,sim) = permute(ival(1,2,:),[3 2 1]);
            end
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            s.processed.splithalf.wholemap.p = 1-nanmean(bsxfun(@gt,...
                permute(s.processed.splithalf.wholemap.val(1,2,:),[3 2 1]),[null(repmat(include,[1 nsims]))]'),2);


            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.wholemap.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/WholeMap/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.wholemap.p <= 0.05) ...
                length(s.processed.splithalf.wholemap.val(1,2,:)) nanmean(s.processed.splithalf.wholemap.p <= 0.05).*100]);
        end

        
        if ismember({'roomxdoor'},varargin)
            didChange = true;
            %%% In-room split half reliability, include doorway

            vel = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
            vel = imfilter(vel,fspecial('gauss',[1 30],10),'same','replicate');
            
            partitions = [1; find(diff(s.processed.validTraceFrames(:,2))>300); length(s.processed.validTraceFrames(:,2))+1];
            
            for pi = 1:length(partitions)-1
                mask = false(1,length(s.processed.p(1,:)));
                mask(partitions(pi):partitions(pi+1)-1) = true;
                allMasks{pi} = [mask & vel>velThresh];
            end
            [a b c ivals] = getMatchedMapsNMasks(s.processed.p,s.processed.trace,allMasks);

            tmp = ivals(1:end/2,end/2+1:end,:);
            tmp = permute(tmp,[3 2 1]);
            allComp = reshape(tmp,length(tmp(:,1,1)),[]);


            allComp = sort(allComp,2,'descend');
            actual = allComp;
            s.processed.splithalf.vals = allComp;


            null = nan(length(s.processed.trace(:,1)),length(allMasks),nsims);
            fprintf('\t\tWithin room by doorways, computing null...  ')
            tic
            P = s.processed.p;
            T = s.processed.trace;
            clear ivals
            parfor sim = 1:nsims
                gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);

                [map samp allComp ivals] = getMatchedMapsNMasks(P,gT,allMasks);

                tmp = ivals(1:end/2,end/2+1:end,:);
                tmp = permute(tmp,[3 2 1]);
                allComp = reshape(tmp,length(tmp(:,1,1)),[]);

                null(:,:,sim) = allComp;
            end
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            null = sort(null,2,'descend');
            s.processed.splithalf.null = null;    
            tmp = nanmax(null,[],2);
            s.processed.splithalf.p = 1-nanmean(bsxfun(@gt,...
                nanmax(s.processed.splithalf.vals,[],2)',permute(tmp,[3 1 2])))';

            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.p <= 0.05) ...
                length(s.processed.splithalf.vals(:,1)) nanmean(s.processed.splithalf.p <= 0.05).*100]);
        end
        
        if didChange
            save(p{1},'-struct','s','-v7.3');
        end
    end
end

function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end

function sh = help_splithalf(p,t)

    v = [0 sqrt(sum(diff(p,[],2).^2))].*30;
    m1 = mkTraceMaps(p,t,1:length(p(1,:))<length(p(1,:))./2);
    m2 = mkTraceMaps(p,t,1:length(p(1,:))>=length(p(1,:))./2);

    rm1 = nan(length(m1(1,:,1)).*length(m1(:,1,1)),length(m1(1,1,:)));
    rm2 = nan(length(m1(1,:,1)).*length(m1(:,1,1)),length(m1(1,1,:)));
    for k = 1:length(m1(1,1,:))
        tmp = m1(:,:,k);
        rm1(:,k) = tmp(:);
        tmp = m2(:,:,k);
        rm2(:,k) = tmp(:);
    end

    isBad = isnan(rm1) | isnan(rm2);

    m1t = reshape(rm1(~isBad),[nansum(~isBad(:,1)) length(m1(1,1,:))]);
    m2t = reshape(rm2(~isBad),[nansum(~isBad(:,1)) length(m1(1,1,:))]);

    xc = (corr(m1t,m2t));
    sh = xc(logical(eye(size(xc))));
end