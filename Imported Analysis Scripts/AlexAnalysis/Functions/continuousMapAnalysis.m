function continuousMapAnalysis(paths,doPlot)
    
    clc
    close all
    drawnow
    
    pause_thresh = 2;
    cell_thresh_for_burst = 0.1;
    doPreplay = false;
    doSubBins = false;
    subBinSize = 60;
    doCOM = false;
    if nargin < 2
        doPlot = false;
    end

    warning off all
     %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    labels = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
    
    envSize = [17 17];
    includeLags = inf;
    pThresh = 0.025;
    envLabel = [{'sq1'} {'sq2'} {'sq3'} ...
        {'g3'} {'g2'} {'g1'}];
    doComps = [{'sq1'} {'sq1'}; {'sq1'} {'g1'}; {'sq1'} {'sq2'}; {'g1'} {'sq2'}; ...
        {'sq1'} {'sq3'}; {'g1'} {'sq3'}; {'sq1'} {'g3'}; {'g1'} {'g3'}; ...
        {'sq1'} {'g2'}; {'g1'} {'g2'}; {'sq1'} {'g1'}; {'g1'} {'g1'}]; %%% Complete curve
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        doAl = help_getAlignmentID(s.alignment,length(isM),paths(isM));
        alignMap = s.alignment(doAl).alignmentMap;
        am = repmat({[]},[1 length(sessions)]);
        amfr = repmat({[]},[1 length(sessions)]);
        isPC = repmat({[]},[1 length(sessions)]);
        COMs = repmat({[]},[1 length(sessions)]);
        aos = [];
        envs = [];
        tic
        fprintf(['\t\tPreloading Data... '])
        for si = 1:length(sessions)
            s = load(sessions{si},'processed','exclude');
            slashInds = find(ismember(sessions{si},'/'));
            gT = s.processed.trace;
            if isfield(s.processed,'exclude')
                gT = gT(s.processed.exclude.SFPs,:);
            end
            amfr{si} = nanmean(gT,2);
            v = [0 sqrt(nansum(diff(s.processed.p,[],2).^2,1))].*30;
            [m os] = mkTraceMaps(s.processed.p,gT,v>=pause_thresh,envSize);
            aos = cat(3,aos,os);
            am{si} = m;
            isPC{si} = s.processed.splithalf.wholemap_unmatched.p <= pThresh;
            if isfield(s.processed,'exclude')
                isPC{si} = isPC{si}(s.processed.exclude.SFPs,:);
            end
            envs = [envs; {lower(sessions{si}(find(ismember(sessions{si},'_'),1,'last')+1:end-4))}];
            if doPreplay
                [preVecs{si} pvals] = getBurstVectors(s.processed.p,gT, ...
                    pause_thresh,cell_thresh_for_burst,isPC{si});
            end
            if doCOM
                COMs{si} = getCOMs(m);
            end
        end  
        durat = toc;
        fprintf([num2str(durat) ' s\n']);
        
        um = nan([envSize length(alignMap{1}(:,1)) length(sessions)]);
        for si = 1:length(sessions)
%             alignMap{1}(~isPC{si},si) = 0; %%% Remove the non-place cells
            um(:,:,alignMap{1}(:,si)~=0,si) = am{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));
            if doPreplay
                tmp = nan(length(alignMap{1}(:,1)),length(preVecs{si}(1,:)));
                tmp(alignMap{1}(:,si)~=0,:) = preVecs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si),:);
                preVecs{si} = tmp;
            end
            if doCOM
                tmp = nan(length(alignMap{1}(:,1)),length(COMs{si}(1,:)));
                tmp(alignMap{1}(:,si)~=0,:) = COMs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si),:);
                COMs{si} = tmp;
            end
        end  
        
        if doPreplay
            tic
            iter = 0;
            strLength = 0;
            fprintf(['\t\tSVM preplay vec classification... '])
            preVecSVMs = nan([length(sessions) length(sessions) length(sessions)]);
            for vecI = 1:length(sessions)
                for si = 1:length(sessions)  
                    for sj = si+1:length(sessions)   
                        iter = iter+1;
                        fprintf(repmat('\b',[1 strLength]));
                        str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2).*length(sessions))]);
                        fprintf(str);
                        strLength = length(str);
                
                        
                        cv1 = permute(um(:,:,:,si),[3 1 2]);
                        cv1 = reshape(cv1,[length(cv1(:,1,1)) prod(size(um(:,:,1,1)))]);
                        cv2 = permute(um(:,:,:,sj),[3 1 2]);
                        cv2 = reshape(cv2,[length(cv2(:,1,1)) prod(size(um(:,:,1,1)))]);
                        goodC = ~all(isnan(cv1),2) & ~all(isnan(cv2),2) & ...
                            ~all(isnan(preVecs{vecI}),2);
                        cv1(~goodC,:) = [];
                        cv2(~goodC,:) = [];
                        goodP = ~all(isnan(cv1),1) & ~all(isnan(cv2),1);
                        cv1(:,~goodP) = [];
                        cv2(:,~goodP) = [];
                        
                        tsvm = svmtrain([cv1 cv2]', ...
                            [ones(1,length(cv1(1,:))) 2.*ones(1,length(cv2(1,:)))]');
                        [predLabel] = svmclassify(tsvm,preVecs{vecI}(goodC,:)');
                        preVecSVMs(si,sj,vecI) = nanmean((predLabel-1));
                    end
                end
            end
            durat = toc;
            fprintf([num2str(durat) ' s\n']);
        end
        
        kDeformFit = nan([length(sessions) length(sessions) length(alignMap{1}(:,1))]);
        kDeformSim = nan([length(sessions) length(sessions) length(alignMap{1}(:,1))]);
        sim = nan([length(sessions) length(sessions) length(alignMap{1}(:,1))]);
        kCrossCOM = repmat({[]},length(sessions));
        iter = 0;
        strLength = 0;
        fprintf(['\n\t\tComputing pairwise map comparisons: '])
        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                iter = iter+1;
                fprintf(repmat('\b',[1 strLength]));
                str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2))]);
                fprintf(str);
                strLength = length(str);
                
                tmp1 = um(:,:,:,si);
                tmp2 = um(:,:,:,sj);
                
% % %                 isOff1 = permute(all(all(isnan(tmp1),1),2),[3 1 2]);
% % %                 isOff2 = permute(all(all(isnan(tmp2),1),2),[3 1 2]);
% % %                 
% % %                 filler1 = nan(size(tmp1(:,:,1)));
% % %                 filler1(~all(isnan(tmp1),3)) = 0;
% % %                 tmp1(:,:,isOff1) = repmat(filler1,[1 1 nansum(isOff1)]);
% % %                 
% % %                 filler2 = nan(size(tmp2(:,:,1)));
% % %                 filler2(~all(isnan(tmp2),3)) = 0;
% % %                 tmp2(:,:,isOff2) = repmat(filler2,[1 1 nansum(isOff2)]);
                
                
                xc = pvxcorr3(tmp1,tmp2,[1 1 1],20);
%                 [a b c] = ind2sub(size(xc),nanmedian(find(xc==nanmax(xc(:)))));
                a = 2;
                b = 2;
                c = 2;
                ivals = xcorr3transform(tmp1,tmp2,[a-ceil(length(xc(:,1,1))./2) ...
                    b-ceil(length(xc(1,:,1))./2) c-ceil(length(xc(1,1,:))./2)]);
                sim(si,sj,:) = ivals;
                
                isOff1 = permute(all(all(isnan(tmp1),1),2),[3 1 2]);
                isOff2 = permute(all(all(isnan(tmp2),1),2),[3 1 2]);

                [kDeformSim(si,sj,~isOff1&~isOff2) kDeformFit(si,sj,~isOff1&~isOff2)] = ...
                    fitDeformation(tmp1(:,:,~isOff1&~isOff2), ...
                    tmp2(:,:,~isOff1&~isOff2)); %x y both neither
                
                %%% Correct for those that are best fit by no change
                [a b] = nanmax([permute(kDeformSim(si,sj,~isOff1&~isOff2),[3 1 2]) ...
                    ivals(~isOff1&~isOff2)],[],2);
                tind = find(~isOff1&~isOff2);
                kDeformSim(si,sj,tind) = a;
                kDeformFit(si,sj,tind(b==2)) = 4;
                
                
%                 figure(1)
%                 scatter(permute(sim(si,sj,~isOff1&~isOff2),[3 1 2]), ...
%                     permute(kDeformSim(si,sj,~isOff1&~isOff2),[3 1 2]))
%                 axis equal
%                 axis square
%                 drawnow
                if doCOM
                    kCrossCOM{si,sj} = [COMs{si} COMs{sj}];
                end
            end
        end
        
        slashInds = find(ismember(paths{1},'/'));
        root = ['Plots/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
        slashInds = find(ismember(upiece{mi},'/'));
        root = [root '/' upiece{mi}(slashInds(end)+1:end) '/ContinuousAnalyses'];
        
        
        percentDeforms = nan([length(sessions) length(sessions) 4]);
        tmp = kDeformFit;
        tmp(kDeformSim<0.7) = nan;
        for i = 1:4
            percentDeforms(:,:,i) = nansum(tmp==i,3)./nansum(~isnan(tmp),3);
        end
        
        percentDeforms = nan([length(sessions) length(sessions) 4]);
        for i = 1:4
            percentDeforms(:,:,i) = nansum(kDeformFit==i,3)./nansum(~isnan(kDeformFit),3);
        end
        
        if doCOM
            [a condGroup] = ismember(envs,envLabel);
            comVecPlot(groupMat(kCrossCOM,condGroup,false),[root '/COM_Vectors']);
        end
        
%         tmp = repmat({[]},[length(envs)]);
%         for i = 1:length(envs)
%             for j = 1:length(envs)
%                 tmp{i,j} = sim(i,j,~isnan(sim(i,j,:)));
%             end
%         end
%         
%         params = transitionPlot(tmp,envs,doComps, ...
%             [root '/ConditionComparison_Correlations_Cellwise']);
        
        mds2D(kDeformSim,envs,envLabel,[root '/MDS_2D']);
        mds3D(kDeformSim,envs,envLabel,[root '/MDS_3D.gif']);
        
        if doPreplay
            figure
            tmp = (preVecSVMs-0.5);
            step = 0.01;
            cm = [[[0:step:1]'; ones(length(0:step:1),1)] ...
                [[0:step:1 1:-step:0]'] ...
                flipud([[0:step:1]'; ones(length(0:step:1),1)])];
            for gi = 1:length(envLabel)
                imagesc(nanmedian(tmp(:,:,ismember(envs,envLabel(gi))),3));
                colormap(cm);
                caxis([-0.5 0.5])
                saveFig(gcf,[root '/BurstVec_SVM_' envLabel{gi}],[{'tiff'} {'pdf'}])
            end
            
            doC = bsxfun(@times,ismember(envs,{'sq1'}),ismember(envs,{'g1'})');
            s1Tog1 = reshape(tmp(logical(repmat(doC, ...
                [1 1 length(sessions)]))),[nansum(doC(:)) length(sessions)]);
            doC = bsxfun(@times,ismember(envs,{'g1'}),ismember(envs,{'sq1'})');
            g1Tosq1 = -reshape(tmp(logical(repmat(doC, ...
                [1 1 length(sessions)]))),[nansum(doC(:)) length(sessions)]);
            comb = [s1Tog1; g1Tosq1];
            [a b] = nanmax(abs(comb),[],1);
            for i = 1:length(sessions)
                maxbias(i) = comb(b(i),i);
            end
            
            transitionPlot(crossSim,envs,doComps, ...
                [root '/ConditionComparison_PV']);
        end
        
        %%%%%%%%%%%%%%%%%% RDM ANALYSIS BEGIN %%%%%%%%%%%%%%
        
        features = [];
        for i = 1:length(envLabel)
            envI = find(ismember(envs,envLabel(i)));
            for k = 1:length(envI)
                features(envI(k),:) = [k i];
            end
        end
        lagMat = abs(bsxfun(@minus,[1:length(sim(:,1,1))],[1:length(sim(:,1,1))]'));
        iterMat = abs(bsxfun(@minus,features(:,1),features(:,1)'));
        envMat = abs(bsxfun(@minus,features(:,2),features(:,2)'));
        attractorMat = abs(bsxfun(@minus,features(:,2)>=4,[features(:,2)>=4]'));
% % %         morphPoint = [4 4 5 5 6 6];
% % %         vec = [];
% % %         for gi = 1:6
% % %             vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
% % %         end
% % %         attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
% % %             [features(:,2)'>=vec(1:length(features(:,2)))]'));
        
%         lagMat = lagMat-nanmean(lagMat(:));
%         iterMat = iterMat-nanmean(iterMat(:));
%         envMat = envMat-nanmean(envMat(:));
        fits = nan(length(sim(1,1,:)),6);
        subFits = nan(length(sim(1,1,:)),3);
        for k = 1:length(sim(1,1,:))
            tmp = kDeformSim(:,:,k);
            if all(isnan(tmp(:)))
                continue
            end
            fits(k,:) = [corr(tmp(~isnan(tmp)),lagMat(~isnan(tmp)),'type','kendall') ...
                corr(tmp(~isnan(tmp)),attractorMat(~isnan(tmp)),'type','kendall') ...
                corr(tmp(~isnan(tmp)),envMat(~isnan(tmp)),'type','kendall') ...
                corr(tmp(~isnan(tmp)),attractorMat(~isnan(tmp)).*lagMat(~isnan(tmp)),'type','kendall') ...
                corr(tmp(~isnan(tmp)),envMat(~isnan(tmp)).*lagMat(~isnan(tmp)),'type','kendall') ...
                nansum(~isnan(tmp(:)))];
            subTmp = tmp-fits(k,2).*attractorMat;
            subFits(k,:) = [corr(subTmp(~isnan(tmp)),iterMat(~isnan(tmp)),'type','kendall') ...
                corr(subTmp(~isnan(tmp)),attractorMat(~isnan(tmp)),'type','kendall') ...
                corr(subTmp(~isnan(tmp)),envMat(~isnan(tmp)),'type','kendall')];
            
        end
        fits = fits(:,[1 2 3 4 6]);

        minSamples = nchoosek(length(sessions),2).*0.5;
        colors = bsxfun(@times,[0 0.25 1],-fits(fits(:,end)>minSamples,2));
        colors = colors+bsxfun(@times,[1 0 0],-fits(fits(:,end)>minSamples,1));
        colors(colors<0) = 0;
        colors = colors.^(1./4);
        colors = [colors./nanmax(colors(:))];
        figure
        set(gcf,'position',[50 50 350 350])
        scatter(-fits(fits(:,end)>minSamples,2), ...
            -fits(fits(:,end)>minSamples,1),...
            30,colors,'filled')
        hold on
        xlabel('Attractor RDM Fit (Kendall''s \tau)')
        ylabel('Time RDM Fit (Kendall''s \tau)')
        set(gca,'xlim',[-0.4 0.8],'ylim',[-0.4 0.8])
        axis square
        plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
        plot([0 0],get(gca,'ylim'),'linestyle','--','color','k')
        saveFig(gcf,[root '/RDM_Analysis_Cellwise'],[{'pdf'} {'tiff'}])
   
        [goodInds] = find(fits(:,end)>minSamples);
        figure
        set(gcf,'position',[50 50 700 700])
        subplot(2,2,1)
        imagesc(-lagMat)
        axis off
        subplot(2,2,2)
        imagesc(-attractorMat)
        axis off
        subplot(2,2,3)
        [a best] = sort(abs(fits(fits(:,end)>minSamples,1)),'descend');
        ex = squarify(kDeformSim(:,:,(goodInds(best(2)))));
        imagesc(ex)
        alpha(double(~isnan(ex)))
        axis off
        subplot(2,2,4)
        [a best] = sort(abs(fits(fits(:,end)>minSamples,2)),'descend');
        ex = squarify(kDeformSim(:,:,(goodInds(best(1)))));
        imagesc(ex)
        alpha(double(~isnan(ex)))
        axis off
        saveFig(gcf,[root '/RDM_Analysis_Cellwise_Examples'],[{'pdf'} {'tiff'}])
        
        %%%%%%%%%%%%%%%%% RDM ANALYSIS END %%%%%%%%%%%%%%%%%
        
        %%%%%%%%% GEO FITS VERSUS DRIFT / ATTRACTOR LOADINGS
        
%         figure
%         set(gcf,'position',[50 50 1400 700])
%         for i =1:4
%             geofits{i} = permute(nansum(nansum(kDeformFit==i,1),2)./ ...
%                 nansum(nansum(~isnan(kDeformFit),1),2),[3 4 1 2]);
%             subplot(2,4,i)
%             scatter(geofits{i}(fits(:,end)>minSamples),-fits(fits(:,end)>minSamples,1))
%             lsline
%             subplot(2,4,i+4)
%             scatter(geofits{i}(fits(:,end)>minSamples),-fits(fits(:,end)>minSamples,2))
%             lsline
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         mds3D(sim,envs,envLabel,[root '/MDS_3D.gif']);
%         mds3D(sim(:,:,alignedSessions<=16),envs,envLabel,[root '/MDS_3D_LessRegistered.gif']);
%         mds3D(sim(:,:,alignedSessions>16),envs,envLabel,[root '/MDS_3D_MoreRegistered.gif']);
%         
%         mkSessionComps(sim,envs,doComps,[root '/ConditionComparison_Cellwise']);
%         step = 5;
%         for i = step:step:35
%             mkSessionComps(sim(:,:,(alignedSessions>(i-step))&(alignedSessions<=(i))),envs,doComps, ...
%                 [root '/ConditionComparison_Cellwise_Registered_' num2str(i-step) '_to_' num2str(i)]);
%         end
%         tmp = repmat({[]},[size(kDeformSim(:,:,1))]);
%         for i = 1:length(kDeformSim(1,:,1))
%             for j = 1:length(kDeformSim(1,:,1))
%                 tmp{i,j} = permute(sim(i,j,~isnan(kDeformSim(i,j,:))),[3 1 2]);
%             end
%         end
% 
%         for gi = 1:6:length(kDeformSim(1,:,1))
%             params = transitionPlot(tmp(gi:nanmin(gi+7,length(kDeformSim(1,:,1))), ...
%                 gi:nanmin(gi+7,length(kDeformSim(1,:,1)))), ...
%                 envs(gi:nanmin(gi+7,length(kDeformSim(1,:,1)))),doComps, ...
%                 [root '/ConditionComparison_Correlations_Cellwise_Sequence' num2str(((gi-1)./6) + 1)]);
%         end
      
    end
end