
function alternatingSessionTests()
    load('combinedData_Exp1');
    
    clc
    
    doColors = [1 1 1; ...
        0.9 0.0 0.9; ...
        0.0 0.9 0.9];
    
    cstack = [ ... %repmat(permute(doColors(1,:),[3 1 2]),[12 36]); ...
        repmat(permute(doColors(2,:),[3 1 2]),[12 36]); ...
        repmat(permute(doColors(3,:),[3 1 2]),[12 36])];
    
    %%% Reliability constaint
    
    for mi = 1:length(dat)
        dat(mi).rely = repmat({[]},[6 6]);
        dat(mi).doskip = false(6,6);
        fprintf(['\n\tMouse:  ' num2str(mi) '\n'])
        for di = 1:length(dat(mi).maps.rooms(1,:))
            for dj = di+1:length(dat(mi).maps.rooms(1,:))
                
                if isempty(dat(mi).registration{di,dj})
                    dat(mi).doskip(di,dj) = true;
                    dat(mi).doskip(di,dj) = true;
                    continue
                end
                
                vals = nan(1,4);
                for rot = [0 2]
                    doM1 = dat(mi).maps.overall{:,di};
                    doM2 = dat(mi).maps.overall{:,dj};
                    
                    doM1 = doM1(:,:,dat(mi).registration{di,dj}(:,1));
                    doM2 = doM2(:,:,dat(mi).registration{di,dj}(:,2));
                    
                    doM2 = rot90(doM2,rot);
                    
                    isGood = ~isnan(doM1)&~isnan(doM2);
                    vals(rot+1) = corr(doM1(isGood),doM2(isGood));
                end
            
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
    
        %%% Sampling
    
    
    figure(1)
    set(gcf,'position',[50 50 250.*length(dat) 125.*length(dat)])
    for mi = 1:length(dat)
        for si = 1:length(dat(mi).path)
            subplot(length(dat),length(dat(mi).path),(mi-1).*length(dat(mi).path)+si)
%             if dat(mi).doskip(si)
%                 plot(dat(mi).path{si}(2,:),dat(mi).path{si}(1,:),'color',[0.6 0.6 0.6],'linewidth',1)
%             else
                plot(dat(mi).path{si}(2,:),dat(mi).path{si}(1,:),'color','k','linewidth',1)
%             end
            set(gca,'ydir','reverse')
            axis equal
            axis off
        end
    end
    root = 'Plots/Summary/';
    saveFig(gcf,[root 'Sampling'],'pdf');
    close all
    drawnow
    
    %% Plot maps combined across days
    
    for mi = 1:length(dat)
        simXLocation = repmat({[]},[6 6]);
        animalSims = [];
        animalSimR = [];
        for di = 1:1:length(dat(mi).maps.rooms(1,:))
            for dj = di+1:length(dat(mi).maps.rooms(1,:))
                if dat(mi).doskip(di,dj)
                    continue
                end
                gT = [dat(mi).trace{di}(dat(mi).registration{di,dj}(:,1),:) ...
                    dat(mi).trace{dj}(dat(mi).registration{di,dj}(:,2),:)];
                gT = gT(dat(mi).rely{di,dj},:);
                doP = [dat(mi).path{di} dat(mi).path{dj}];
                m = mkTraceMaps([dat(mi).path{di} dat(mi).path{dj}],gT,[]);

                doK = [9 8];
                for part = 0:floor(length(gT(:,1))/prod(doK))

                    step = 30;

                    close all

                    com = [dat(mi).mostrecent.common{di} dat(mi).mostrecent.common{dj}];
                    uni = [dat(mi).mostrecent.unique{di} dat(mi).mostrecent.unique{dj}];
                    isInRoom = [dat(mi).isinroom{di} dat(mi).isinroom{dj}];

                    clear maps
                    for room = 1:3
                        maps.rooms{room} = mkTraceMaps( ...
                            doP(:,[isInRoom(room,:)]),...
                            gT(:,[isInRoom(room,:)]),[],[12 12]);

                        maps.common{room} = ...
                            mkTraceMaps(doP(:,[isInRoom(room,:)]),...
                            gT(:,[isInRoom(room,:)]),com(isInRoom(room,:)),[12 12]);

                        maps.unique{room} = ...
                            mkTraceMaps(doP(:,[isInRoom(room,:)]),...
                            gT(:,[isInRoom(room,:)]),uni(isInRoom(room,:)),[12 12]);
                    end

                    m = nan(36,36,length(maps.rooms{1}(1,1,:)));
                    for r = 1:3
                        m(1:12,(r-1).*12+1:(r).*12,:) = maps.rooms{r};
                        m([13:24],(r-1).*12+1:(r).*12,:) = maps.common{r};
                        m([25:36],(r-1).*12+1:(r).*12,:) = maps.unique{r};
                    end

                    figure(1)
                    set(gcf,'position',[50 50 750 800])
                    for k = 1:prod(doK)
                        if part.*prod(doK)+k > length(gT(:,1))
                            break
                        end
                        subplot(doK(1),doK(2),k)
                        imagesc(m(1:12,:,part.*prod(doK)+k))
                        colormap('jet')
                        alpha(double(~isnan(m(1:12,:,part.*prod(doK)+k))))

                        axis equal
                        axis off
                    end

                    outP = ['Plots/TracePlots_Combined/ThreeRooms/Animal_' num2str(mi)  '/Day_' num2str(di) '_Partition_' num2str(part+1)];
                    saveFig(gcf,outP,'tiff')
                    saveFig(gcf,outP,'pdf')
                    close all
                    drawnow

                    close all


    %                 m = m./repmat(nanmax(nanmax(m,[],2),[],1),[size(m(:,:,1))]);

                    figure(1)
                    set(gcf,'position',[50 50 750 800])
                    for k = 1:prod(doK)
                        if part.*prod(doK)+k > length(gT(:,1))
                            break
                        end
                        subplot(doK(1),doK(2),k)
                        imagesc(m(13:end,:,part.*prod(doK)+k));
                        colormap('jet')
                        alpha(double(~isnan(m(13:end,:,part.*prod(doK)+k))))

                        axis equal
                        axis off
                    end

                    outP = ['Plots/TracePlots_Combined_Differentiated/ThreeRooms/Animal_' num2str(mi)  '/Day_' num2str(di) '_' num2str(dj) '_Partition_' num2str(part+1)];
                    saveFig(gcf,outP,'tiff')
                    saveFig(gcf,outP,'pdf')
                    close all
                    drawnow
                end                
            end
        end
    end
    
    %% Room map similarity, ignore doors

%     allSim = [];
%     for mi = 1:length(dat)
%         for di = 1:2:length(dat(mi).maps.rooms(1,:))
%             doM1 = cat(4,dat(mi).maps.rooms{:,di});
%             doM2 = cat(4,dat(mi).maps.rooms{:,di+1});
%             doM1 = doM1(:,:,dat(mi).isAligned{di},:);
%             doM2 = doM2(:,:,dat(mi).isAligned{di+1},:);
%             
%             sim = nan(3,3);
%             for i = 1:length(doM1(1,1,1,:))
%                 for j = 1:length(doM2(1,1,1,:))
%                     m1 = doM1(:,:,:,i);
%                     m2 = doM2(:,:,:,j);
%                     isGood = ~isnan(m1)&~isnan(m2);
%                     sim(i,j) = corr(m1(isGood),m2(isGood));
%                 end
%             end
%             allSim = cat(3,allSim,sim);
%         end
%     end
    
    %%% Room map similarity, divide by doors
    
    allSim = [];
    allSim2 = [];
    allSimR = [];
    simXDay = [{[]} {[]} {[]}];
    simXPixel = [{[]} {[]} {[]}];
    nullAllSim = [];
    for mi = 1:length(dat)
        simXLocation = repmat({[]},[6 6]);
        animalSims = [];
        animalSimR = [];
        nullAnimalSims = [];
        for di = 1:length(dat(mi).maps.rooms(1,:))
            for dj = di+1:length(dat(mi).maps.rooms(1,:))
                
                if dat(mi).doskip(di,dj)
                    continue
                end

                doM1 = cat(4,dat(mi).maps.common{:,di});
                doM2 = cat(4,dat(mi).maps.common{:,dj});
                doM1 = cat(4,doM1,dat(mi).maps.unique{:,di});
                doM2 = cat(4,doM2,dat(mi).maps.unique{:,dj});

                doM1 = doM1(:,:,dat(mi).registration{di,dj}(:,1),:);
                doM2 = doM2(:,:,dat(mi).registration{di,dj}(:,2),:);

                doM1 = doM1(:,:,dat(mi).rely{di,dj},:);
                doM2 = doM2(:,:,dat(mi).rely{di,dj},:);

%                 scalar1 = repmat(nanmax(nanmax(doM1,[],2),[],1),[size(doM1(:,:,1))]);
%                 scalar1(scalar1==0) = 1; %%% Correct so don't try to divide by zero
%                 scalar2 = repmat(nanmax(nanmax(doM2,[],2),[],1),[size(doM2(:,:,1))]);
%                 scalar2(scalar2==0) = 1; %%% Correct so don't try to divide by zero

    %             doM1 = doM1./scalar1;
    %             doM2 = doM2./scalar2;
                nsims = 100;
                sim = nan(6,6);
                null = nan(6,6,nsims);
                simR = nan(6,6,4);
                for i = 1:length(doM1(1,1,1,:))
                    for j = 1:length(doM2(1,1,1,:))
                        m1 = doM1(:,:,:,i);
                        m2 = doM2(:,:,:,j);
                        isGood = ~isnan(m1)&~isnan(m2);
                        if ~any(isGood(:))
                            continue
                        end
                        sim(i,j) = corr(m1(isGood),m2(isGood));

                        simX = nan(size(m1(:,:,1)));
                        for qi = 1:length(m1(:,1,1))
                            for qj = 1:length(m1(1,:,1))
                                simX(qi,qj) = corr(permute(m1(qi,qj,:),[3 2 1]),...
                                    permute(m2(qi,qj,:),[3 2 1]));
                            end                        
                        end
                        simXLocation{i,j} = cat(3,simXLocation{i,j},simX);

                        for r = 0:3
                            tmpM2 = imrotate(m2,r.*90);
                            isGood = ~isnan(m1)&~isnan(tmpM2);
                            if ~any(isGood(:))
                                continue
                            end
                            simR(i,j,r+1) = corr(m1(isGood),tmpM2(isGood));
                        end

    %                     for si = 1:nsims
    %                         m1 = doM1(:,:,:,i);
    %                         m2 = doM2(:,:,:,j);
    %                         isGood = ~isnan(m1)&~isnan(m2);
    %                         if ~any(isGood(:))
    %                             continue
    %                         end
    %                         m1 = m1(:,:,randperm(length(m1(1,1,:))));
    %                         null(i,j,si) = corr(m1(isGood),m2(isGood));
    %                     end
                    end
                end
%                 allSim2 = cat(3,allSim2,sim);
                animalSims = cat(3,animalSims,sim);
%                 animalSimR = cat(4,animalSimR,simR);
%                 simXDay{(di-1)./2 + 1} = cat(3,simXDay{(di-1)./2 + 1},sim);
            end
        end
        simXLocation = cellfun(@nanmean,simXLocation,...
            num2cell(ones(size(simXLocation)).*3),'uniformoutput',false);

        simXPixel{1} = cat(3,simXPixel{1},nanmean(cat(3,simXLocation{1,4},simXLocation{4,1}),3));
        simXPixel{2} = cat(3,simXPixel{2},nanmean(cat(3,simXLocation{2,5},simXLocation{5,2}),3));
        simXPixel{3} = cat(3,simXPixel{3},nanmean(cat(3,simXLocation{3,6},simXLocation{6,3}),3));

        allSim = cat(3,allSim,nanmean(animalSims,3));        
        allSimR = cat(4,allSimR,nanmean(animalSimR,4)); 
    end
    
    
% %     %%%% Remapping quant for active place cells
% %     
% %     allSimXCell = [];
% %     allSimRXCell = [];
% %     for mi = 1:length(dat)
% %         animalSims = [];
% %         animalSimsR = [];
% %         for di = 1:2:length(dat(mi).maps.rooms(1,:))
% %             if dat(mi).doskip(di)
% %                 continue
% %             end
% %             
% %             doM1 = cat(4,dat(mi).maps.common{:,di});
% %             doM2 = cat(4,dat(mi).maps.common{:,di+1});
% %             doM1 = cat(4,doM1,dat(mi).maps.unique{:,di});
% %             doM2 = cat(4,doM2,dat(mi).maps.unique{:,di+1});
% %             
% %             doM1 = doM1(:,:,dat(mi).isAligned{di},:);
% %             doM2 = doM2(:,:,dat(mi).isAligned{di+1},:);
% % 
% %             doM1 = doM1(:,:,dat(mi).rely{di},:);
% %             doM2 = doM2(:,:,dat(mi).rely{di+1},:);
% %             
% % %             length(doM1(1,1,:,1))
% %             simR = nan(6,6,length(doM1(1,1,:,1)),4);
% %             sim = nan(6,6,length(doM1(1,1,:,1)));
% %             for i = 1:length(doM1(1,1,1,:))
% %                 for j = 1:length(doM2(1,1,1,:))
% %                     for k = 1:length(doM1(1,1,:,1))
% %                         m1 = doM1(:,:,k,i);
% %                         m2 = doM2(:,:,k,j);
% %                         isGood = ~isnan(m1)&~isnan(m2);
% %                         if ~any(isGood(:))
% %                             continue
% %                         end
% %                         try
% %                             sim(i,j,k) = corr(m1(isGood),m2(isGood));
% %                         end
% %                         
% % %                         for r = 0:3
% % %                             tmpM2 = imrotate(m2,r.*90);
% % %                             isGood = ~isnan(m1)&~isnan(tmpM2);
% % %                             if ~any(isGood(:))
% % %                                 continue
% % %                             end
% % %                             simR(i,j,k,r+1) = corr(m1(isGood),tmpM2(isGood));
% % %                         end
% %                     end
% %                 end
% %             end
% %             animalSims = cat(3,animalSims,(sim));
% % %             animalSimsR = cat(3,animalSimsR,simR);            
% %             simXDay{(di-1)./2 + 1} = cat(3,simXDay{(di-1)./2 + 1},sim);
% %         end
% %         
% %         allSimXCell = cat(3,allSimXCell,animalSims);   
% % %         allSimRXCell = cat(3,allSimRXCell,animalSimsR);
% %     end
% % %     allSim = allSimXCell;
% % 
% % %     [a b] = nanmax(allSimRXCell,[],4);
% % %     hist(permute(cat(3,b(1,4,:)==2,b(4,1,:)),[3 1 2]))
% % %     hist(permute(cat(3,b(2,5,:),b(5,2,:)),[3 1 2]))
% % %     hist(permute(cat(3,b(3,6,:),b(6,3,:)),[3 1 2]))
    
    %% Make plots

    close all
    figure(1)
    set(gcf,'position',[50 50 425 300])
    imagesc(nanmean(allSim(:,:,1:end),3))
%     colormap parula
%     colormap(circshift([linspace(0,1,256)' ...
%         [linspace(0,1,128) linspace(1,0,128)]' ...
%         [linspace(0,1,128) linspace(1,0,128)]'],[0 3]))
%     colormap(circshift([linspace(0,1,256)' ...
%         [linspace(0,0.35,128) linspace(0.35,0,128)]' ...
%         [linspace(0,0.35,128) linspace(0.35,0,128)]'],[0 0]))
    colormap(circshift([linspace(0,1,256)' ...
        [linspace(0,1,128) ones(1,128)]' ...
        linspace(0,1,256)'],[0 -2]))
    caxis([-0.0 0.5])
    colorbar
    axis equal
    axis off
    figure(2)
    set(gcf,'position',[400 50 125 200])
    mask = false(6,6);
    mask(1:7:15) = true;
    v = nan(length(allSim(1,1,:)),3);
    v(:,1) = help_getMaskedVals(allSim,mask);
    v(:,2) = help_getMaskedVals(allSim,circshift(mask,[3 3]));
    v(:,3) = help_getMaskedVals(allSim,circshift(mask,[0 3])|circshift(mask,[3 0]));
    av{1} = v;
    mkGraph(v)
    figure(3)
    set(gcf,'position',[550 50 125 200])
    mask = false(6,6);
    mask([2 3 9 7 13 14]) = true;
    v = nan(length(allSim(1,1,:)),2);
    v(:,1) = help_getMaskedVals(allSim,mask);
    v(:,2) = help_getMaskedVals(allSim,circshift(mask,[3 3]));
    mkGraph(v)
    av{2} = v;
    figure(4)
    set(gcf,'position',[700 50 125 200])
    mask = false(6,6);
    mask([4 19]) = true;
    v = nan(length(allSim(1,1,:)),3);
    v(:,1) = help_getMaskedVals(allSim,mask);
    v(:,2) = help_getMaskedVals(allSim,circshift(mask,[1 1]));
    v(:,3) = help_getMaskedVals(allSim,circshift(mask,[2 2]));
    mkGraph(v) 
    av{3} = v;
    
    figure(5)
    set(gcf,'position',[900 50 250 200])
    v = nan(length(allSim(1,1,:)),4);
    mask = false(6,6);
    mask([2 9 7 14]) = true;
    v(:,1) = help_getMaskedVals(allSim,mask);
    v(:,3) = help_getMaskedVals(allSim,circshift(mask,[3 3]));
    mask = false(6,6);
    mask([3 13]) = true;
    v(:,2) = help_getMaskedVals(allSim,mask);
    v(:,4) = help_getMaskedVals(allSim,circshift(mask,[3 3]));
    subplot(1,2,1)
    mkGraph(v(:,1:2))
    set(gca,'ylim',[-0.0 0.25])
    subplot(1,2,2)
    mkGraph(v(:,3:4))
    set(gca,'ylim',[-0.0 0.25])
    av{4} = v;
    
%     figure(6)
%     set(gcf,'position',[1200 50 250 200])
%     v = [];
%     for di = 1:3
%         mask = false(6,6);
%         mask(1:7:15) = true;
%         v(:,di) = help_getMaskedVals(simXDay{di},circshift(mask,[0 3])|circshift(mask,[3 0]));
%     end
%     mkGraph(v)
    
    figure(7)
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
    
    figure(8)
%     set(gcf,'position',[50 450 300 300])
%     doPlot = repmat({[]},[4 3]);
%     for r = 1:4
%         for i = 1:3
%         doPlot{r,i} = nanmean([permute(allSimR(i,i+3,r,:),[4 1 2 3]) ...
%             permute(allSimR(i+3,i,r,:),[4 1 2 3])],2);
%         end
%     end
%     mkGraph(doPlot)
    
    figure(9)
    set(gcf,'position',[450 150 300 150])
    tmp = [];
    for i = 1:4
        tmp = [tmp (dat(i).entrances.unique)-0.5];
    end
    hist(tmp)
    set(gca,'xlim',[-0.25 0.25])
    
    figure(10)
%     set(gcf,'position',[1000 150 250 250])
%     mask = false(6,6);
%     mask(1:7:15) = true;
%     v = [];
%     v(:,1) = help_getMaskedVals(allSimXCell,mask);
%     v(:,2) = help_getMaskedVals(allSimXCell,circshift(mask,[3 3]));
%     v(:,3) = help_getMaskedVals(allSimXCell,circshift(mask,[0 3])|circshift(mask,[3 0]));
%     cumHist(v,[-1:0.025:1])
%     axis square
    
    root = 'Plots/Summary/ThreeRooms/';
    figure(1)
    saveFig(gcf,[root 'RDM'],'pdf');
    figure(2)
    saveFig(gcf,[root 'Reliability'],'pdf');
    figure(3)
    saveFig(gcf,[root 'Repetition'],'pdf');
    figure(4)
    saveFig(gcf,[root 'WithinRoomSimilarity'],'pdf');
    figure(5)
    saveFig(gcf,[root 'RepetitionXDistance'],'pdf');
%     figure(6)
%     saveFig(gcf,[root 'AcrossReliabilityXTime'],'pdf');
    figure(7)
    saveFig(gcf,[root 'SimilarityXPixel'],'pdf');
    figure(8)
    saveFig(gcf,[root 'BestMatchRotation'],'pdf');
    figure(9)
    saveFig(gcf,[root 'DoorwayBehavioralBias'],'pdf');
    figure(10)
    saveFig(gcf,[root 'CumulativeCellDistribution'],'pdf');
    
    %% ALL STATS WRITEOUT
    
    outP = ['Stats_ThreeRooms.txt'];
    fid = fopen(outP,'w');
    [h p ci tstat] = ttest(av{1}(:,1),av{1}(:,2));
    
    fprintf(fid,'\t\t\tSIMILARITY\n');
    fprintf(fid,'\nSimilarity Within-Common vs. Within-Unique:  ');
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    [h p ci tstat] = ttest(av{1}(:,1),av{1}(:,3));
    fprintf(fid,'\nSimilarity Within-Common vs. Across:  ');
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    [h p ci tstat] = ttest(av{1}(:,2),av{1}(:,3));
    fprintf(fid,'\nSimilarity Within-Unique vs. Across:  ');
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    
    [h p ci tstat] = ttest(av{2}(:,1),av{2}(:,2));
    fprintf(fid,'\n\t\t\tREPETITION\n');
    fprintf(fid,'\nCommon vs. Unique:  ');
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    fprintf(fid,'\n\n\t\tACROSS DOOR SIMILARITY BY ROOM\n');
    
    [h p ci tstat] = ttest(av{3}(:,1),av{3}(:,2));
    fprintf(fid,'\nLeft vs. middle:  ');
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    [h p ci tstat] = ttest(av{3}(:,2),av{3}(:,3));
    fprintf(fid,'\nRight vs. middle:  ');
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    [h p ci tstat] = ttest(av{3}(:,1),av{3}(:,3));
    fprintf(fid,'\nLeft vs. Right:  ');
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    fprintf(fid,'\n\n\t\tSIMILARITY ACROSS ROOMS BY DISTANCE\n');
    
    [h p ci tstat] = ttest(av{4}(:,1),av{4}(:,2));
    fprintf(fid,'\nNear vs. Far (Common):  ');
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    [h p ci tstat] = ttest(av{4}(:,3),av{4}(:,4));
    fprintf(fid,'\nNear vs. Far (Unique):  ');
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    fprintf(fid,'\n\n\t\tBEHAVIORAL DOOR PREFERENCE\n');
    [h p ci tstat] = ttest(tmp);
    fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.3f ',tstat.df,tstat.tstat,p);
    
    fclose(fid);
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end















