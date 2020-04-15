
function pairwiseSessionTests_Exp3()
    
    clc
    fprintf('\n')
    load('combinedData_Exp3');
     
%     dat(3) = [];
    %%% Reliability constaint

    for mi = 1:length(dat)
        dat(mi).rely = repmat({[]},[6 6]);
        dat(mi).doskip = false(6,6);
        fprintf(['\n\tMouse:  ' num2str(mi) '\n'])
        
        for di = 1:1:length(dat(mi).maps.rooms(1,:))
            for dj = di+1:1:length(dat(mi).maps.rooms(1,:))
            
                %%% Test all rotations
                vals = [1 0 0 0]; %nan(1,4);
%                 for rot = 0:3
%                     doM1 = dat(mi).maps.overall{:,di};
%                     doM2 = dat(mi).maps.overall{:,dj};
%                     
%                     doM1 = doM1(:,:,dat(mi).registration{di,dj}(:,1));
%                     doM2 = doM2(:,:,dat(mi).registration{di,dj}(:,2));
%                     
%                     doM2 = rot90(doM2,rot);
% 
%                     isGood = ~isnan(doM1)&~isnan(doM2);
%                     vals(rot+1) = corr(doM1(isGood),doM2(isGood));
%                 end

                [blah bestRot] = nanmax(vals);

                doM1 = dat(mi).maps.overall{:,di};
                doM2 = dat(mi).maps.overall{:,dj};

                doM1 = doM1(:,:,dat(mi).registration{di,dj}(:,1));
                doM2 = doM2(:,:,dat(mi).registration{di,dj}(:,2));
                    
                doM2 = rot90(doM2,bestRot-1);
                reliability = nan(length(doM1(1,1,:)),1);
                for k = 1:length(doM1(1,1,:))
                    m1 = doM1(:,:,k);
                    m2 = doM2(:,:,k);
                    isGood = ~isnan(m1)&~isnan(m2);
                    reliability(k) = corr(m1(isGood),m2(isGood));
                end

                nsims = 250;
                null = nan(length(reliability),nsims);
                for si = 1:nsims
                    doM1 = doM1(:,:,randperm(length(doM1(1,1,:))));
                    for k = 1:length(doM1(1,1,:))
                        m1 = doM1(:,:,k);
                        m2 = doM2(:,:,k);
                        isGood = ~isnan(m1)&~isnan(m2);
                        null(k,si) = corr(m1(isGood),m2(isGood));
                    end
                end

                [a b] = sort(null(:));
                thresh = a(round(length(a).*0.95));
                dat(mi).rely{di,dj} = reliability >= thresh;
                dat(mi).rely{di,dj} = reliability >= thresh;
                dat(mi).rotation{di,dj} = bestRot-1;

                fprintf(['\t\t( ' num2str(di) ', ' num2str(dj) '; rot:  ' num2str((bestRot-1).*90) ') Count:  ' num2str(nansum(reliability >= thresh)) ...
                    '\tProportion:  ' num2str(nanmean(reliability >= thresh)) '\n'])

                if nansum(reliability >= thresh) < 20 %%% Toss sessions with fewer than 30 reliable place cells
                    dat(mi).doskip(di,dj) = true;
                    dat(mi).doskip(di,dj) = true;
                else
                    dat(mi).doskip(di,dj) = false;
                    dat(mi).doskip(di,dj) = false;
                end
            end
        end
    end
    
    
% % %     allSim = [];
% % %     allSim2 = [];
% % %     allSimR = [];
% % %     simXDay = [{[]} {[]} {[]}];
% % %     simXPixel = [{[]} {[]} {[]}];
% % %     for mi = 1:length(dat)
% % %         simXLocation = repmat({[]},[6 6]);
% % %         animalSims = [];
% % %         animalSimR = [];
% % %         for di = 1:length(dat(mi).maps.rooms(1,:))
% % %             for dj = di+1:length(dat(mi).maps.rooms(1,:))
% % %                 if dat(mi).doskip(di,dj)
% % %                     continue
% % %                 end
% % %                 
% % %                 sim = nan(3,3);
% % %                 simR = nan(3,3,4);
% % %                 for i = 1:length(dat(mi).door{dj}(:,1))
% % %                     for j = 1:length(dat(mi).door{dj}(:,1))
% % %                 
% % %                         [m1 s1] = mkTraceMaps(dat(mi).path{di}(:,dat(mi).isInRoom{di}),...
% % %                             dat(mi).trace{di}(dat(mi).isInRoom{di}),...
% % %                             dat(mi).door{di}(i,dat(mi).isInRoom{di}),[18 18]);
% % %                         [m2 s2] = mkTraceMaps(dat(mi).path{dj}(:,dat(mi).isInRoom{dj}),...
% % %                             dat(mi).trace{dj}(dat(mi).isInRoom{dj}),...
% % %                             dat(mi).door{dj}(j,dat(mi).isInRoom{dj}),[18 18]);
% % %                         
% % %                         matchSamp = nanmin(s1,s2);
% % %                         
% % %                         [doM1 s1] = mkTraceMaps(dat(mi).path{di}(:,dat(mi).isInRoom{di}),...
% % %                             dat(mi).trace{di}(:,dat(mi).isInRoom{di}),...
% % %                             dat(mi).door{di}(i,dat(mi).isInRoom{di}),[18 18],matchSamp);
% % %                         [doM2 s2] = mkTraceMaps(dat(mi).path{dj}(:,dat(mi).isInRoom{dj}),...
% % %                             dat(mi).trace{dj}(:,dat(mi).isInRoom{dj}),...
% % %                             dat(mi).door{dj}(j,dat(mi).isInRoom{dj}),[18 18],matchSamp);
% % %                         
% % % 
% % %                         doM1 = doM1(:,:,dat(mi).registration{di,dj}(:,1));
% % %                         doM2 = doM2(:,:,dat(mi).registration{di,dj}(:,2));
% % % 
% % %                         doM1 = doM1(:,:,dat(mi).rely{di,dj});
% % %                         doM2 = doM2(:,:,dat(mi).rely{di,dj});
% % % 
% % %                         doM2 = rot90(doM2,dat(mi).rotation{di,dj});
% % % 
% % %                         m1 = doM1;
% % %                         m2 = doM2;
% % %                         isGood = ~isnan(m1)&~isnan(m2);
% % %                         if ~any(isGood(:))
% % %                             continue
% % %                         end
% % %                         sim(i,j) = corr(m1(isGood),m2(isGood));
% % % 
% % %                         for r = 0:3
% % %                             tmpM2 = imrotate(m2,r.*90);
% % %                             isGood = ~isnan(m1)&~isnan(tmpM2);
% % %                             if ~any(isGood(:))
% % %                                 continue
% % %                             end
% % %                             simR(i,j,r+1) = corr(m1(isGood),tmpM2(isGood));
% % %                         end
% % %                     end
% % %                 end
% % %                 allSim2 = cat(3,allSim2,sim);
% % %                 animalSims = cat(3,animalSims,sim);
% % %                 animalSimR = cat(4,animalSimR,simR);
% % %             end
% % %         end
% % %         
% % % %         [a b] = nanmax(simR,[],3);
% % % %         allSim = cat(3,allSim,animalSims);    
% % %         allSim = cat(3,allSim,nanmean(animalSims,3));        
% % % %         allSimR = cat(4,allSimR,nanmean(animalSimR,4)); 
% % %     end

    
    %%% Room map similarity, divide by doors
    
    allSim = [];
    allSim2 = [];
    allSimR = [];
    simXDay = [{[]} {[]} {[]}];
    simXPixel = [{[]} {[]} {[]}];
    for mi = 1:length(dat)
        simXLocation = repmat({[]},[6 6]);
        animalSims = [];
        animalSimR = [];
        for di = 1:length(dat(mi).maps.rooms(1,:))
            for dj = di+1:length(dat(mi).maps.rooms(1,:))
                if dat(mi).doskip(di,dj)
                    continue
                end

                doM1 = cat(4,dat(mi).maps.unique{:,di});
                doM2 = cat(4,dat(mi).maps.unique{:,dj});

                doM1 = doM1(:,:,dat(mi).registration{di,dj}(:,1),:);
                doM2 = doM2(:,:,dat(mi).registration{di,dj}(:,2),:);

                doM1 = doM1(:,:,dat(mi).rely{di,dj},:);
                doM2 = doM2(:,:,dat(mi).rely{di,dj},:);

                doM2 = rot90(doM2,dat(mi).rotation{di,dj});

                sim = nan(3,3);
                simR = nan(3,3,4);
                for i = 1:length(doM1(1,1,1,:))
                    for j = 1:length(doM2(1,1,1,:))
                        m1 = doM1(:,:,:,i);
                        m2 = doM2(:,:,:,j);
                        isGood = ~isnan(m1)&~isnan(m2);
                        if ~any(isGood(:))
                            continue
                        end
                        sim(i,j) = corr(m1(isGood),m2(isGood));

                        for r = 0:3
                            tmpM2 = imrotate(m2,r.*90);
                            isGood = ~isnan(m1)&~isnan(tmpM2);
                            if ~any(isGood(:))
                                continue
                            end
                            simR(i,j,r+1) = corr(m1(isGood),tmpM2(isGood));
                        end
                    end
                end
                allSim2 = cat(3,allSim2,sim);
                animalSims = cat(3,animalSims,sim);
                animalSimR = cat(4,animalSimR,simR);
            end
        end
        
%         [a b] = nanmax(simR,[],3);
        
%         allSim = cat(3,allSim,animalSims);    
        allSim = cat(3,allSim,nanmean(animalSims,3));        
%         allSimR = cat(4,allSimR,nanmean(animalSimR,4)); 
    end

    %% Make plots

    close all
    figure(1)
    set(gcf,'position',[50 50 425 300])
    imagesc(nanmean(allSim(:,:,1:end),3))
%     colormap jet
    colormap(circshift([linspace(0,1,256)' ...
        [linspace(0,1,128) ones(1,128)]' ...
        linspace(0,1,256)'],[0 -2]))
    caxis([0.0 0.5])
    colorbar
    axis equal
    axis off
    
    v = nan(length(allSim(1,1,:)),4);
    mask = false(3,3);
    mask(1:4:end) = true;
    v(:,1) = help_getMaskedVals(allSim,mask);
    
    mask = false(3,3);
    mask(1,2) = true;
    mask(2,1) = true;
    v(:,2) = help_getMaskedVals(allSim,mask);
    
    mask = false(3,3);
    mask(1,3) = true;
    mask(3,1) = true;
    v(:,3) = help_getMaskedVals(allSim,mask);
    
    mask = false(3,3);
    mask(2,3) = true;
    mask(3,2) = true;
    v(:,4) = help_getMaskedVals(allSim,mask);
    
    figure(2)
    set(gcf,'position',[50 450 250 200])
    mkGraph(v)
    set(gca,'ylim',[0 0.5])
    
    
%     av{1} = v;
%     mkGraph(v)
%     figure(3)
%     set(gcf,'position',[550 50 125 200])
%     mask = false(6,6);
%     mask([2 3 9 7 13 14]) = true;
%     v = nan(length(allSim(1,1,:)),2);
%     v(:,1) = help_getMaskedVals(allSim,mask);
%     v(:,2) = help_getMaskedVals(allSim,circshift(mask,[3 3]));
%     mkGraph(v)
%     av{2} = v;
%     figure(4)
%     set(gcf,'position',[700 50 125 200])
%     mask = false(6,6);
%     mask([4 19]) = true;
%     v = nan(length(allSim(1,1,:)),3);
%     v(:,1) = help_getMaskedVals(allSim,mask);
%     v(:,2) = help_getMaskedVals(allSim,circshift(mask,[1 1]));
%     v(:,3) = help_getMaskedVals(allSim,circshift(mask,[2 2]));
%     mkGraph(v) 
%     av{3} = v;
%     
%     figure(5)
%     set(gcf,'position',[900 50 250 200])
%     v = nan(length(allSim(1,1,:)),4);
%     mask = false(6,6);
%     mask([2 9 7 14]) = true;
%     v(:,1) = help_getMaskedVals(allSim,mask);
%     v(:,3) = help_getMaskedVals(allSim,circshift(mask,[3 3]));
%     mask = false(6,6);
%     mask([3 13]) = true;
%     v(:,2) = help_getMaskedVals(allSim,mask);
%     v(:,4) = help_getMaskedVals(allSim,circshift(mask,[3 3]));
%     subplot(1,2,1)
%     mkGraph(v(:,1:2))
%     set(gca,'ylim',[-0.05 0.2])
%     subplot(1,2,2)
%     mkGraph(v(:,3:4))
%     set(gca,'ylim',[-0.05 0.2])
%     av{4} = v;
%     
% %     figure(6)
% %     set(gcf,'position',[1200 50 250 200])
% %     v = [];
% %     for di = 1:3
% %         mask = false(6,6);
% %         mask(1:7:15) = true;
% %         v(:,di) = help_getMaskedVals(simXDay{di},circshift(mask,[0 3])|circshift(mask,[3 0]));
% %     end
% %     mkGraph(v)
%     
%     figure(7)
%     set(gcf,'position',[400 350 900 300])
%     for i = 1:3
%         subplot(1,3,i)
%         h = imagesc(nanmean(simXPixel{i},3));
%         colormap jet
%         caxis([0 0.5])
%         set(h,'alphadata',~isnan(nanmean(simXPixel{i},3)))
%         axis equal
%         axis off
%     end
%     
%     figure(8)
%     set(gcf,'position',[50 450 300 300])
%     doPlot = repmat({[]},[4 3]);
%     for r = 1:4
%         for i = 1:3
%         doPlot{r,i} = nanmean([permute(allSimR(i,i+3,r,:),[4 1 2 3]) ...
%             permute(allSimR(i+3,i,r,:),[4 1 2 3])],2);
%         end
%     end
%     mkGraph(doPlot)
%     
    root = 'Plots/Summary/ShortLong/';
    figure(1)
    saveFig(gcf,[root 'RDM'],'pdf');
    figure(2)
    saveFig(gcf,[root 'Reliability'],'pdf');
%     figure(3)
%     saveFig(gcf,[root 'Repetition'],'pdf');
%     figure(4)
%     saveFig(gcf,[root 'WithinRoomSimilarity'],'pdf');
%     figure(5)
%     saveFig(gcf,[root 'RepetitionXDistance'],'pdf');
% %     figure(6)
% %     saveFig(gcf,[root 'AcrossReliabilityXTime'],'pdf');
%     figure(7)
%     saveFig(gcf,[root 'SimilarityXPixel'],'pdf');
%     figure(8)
%     saveFig(gcf,[root 'BestMatchRotation'],'pdf');
%     
%     %% ALL STATS WRITEOUT
%     
%     outP = ['Stats.txt'];
%     fid = fopen(outP,'w');
%     [h p ci tstat] = ttest(av{1}(:,1),av{1}(:,2));
%     
%     fprintf(fid,'\t\t\tSIMILARITY\n');
%     fprintf(fid,'\nSimilarity Within-Common vs. Within-Unique:  ');
%     fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
%     
%     [h p ci tstat] = ttest(av{1}(:,1),av{1}(:,3));
%     fprintf(fid,'\nSimilarity Within-Common vs. Across:  ');
%     fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
%     
%     [h p ci tstat] = ttest(av{1}(:,2),av{1}(:,3));
%     fprintf(fid,'\nSimilarity Within-Unique vs. Across:  ');
%     fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
%     
%     
%     [h p ci tstat] = ttest(av{2}(:,1),av{2}(:,2));
%     fprintf(fid,'\n\t\t\tREPETITION\n');
%     fprintf(fid,'\nCommon vs. Unique:  ');
%     fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
%     
%     fprintf(fid,'\n\t\tACROSS DOOR SIMILARITY BY ROOM\n');
%     
%     [h p ci tstat] = ttest(av{3}(:,1),av{3}(:,2));
%     fprintf(fid,'\nLeft vs. middle:  ');
%     fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
%     
%     [h p ci tstat] = ttest(av{3}(:,2),av{3}(:,3));
%     fprintf(fid,'\nRight vs. middle:  ');
%     fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
%     
%     [h p ci tstat] = ttest(av{3}(:,1),av{3}(:,3));
%     fprintf(fid,'\nLeft vs. Right:  ');
%     fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
%     
%     fclose(fid);
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end















