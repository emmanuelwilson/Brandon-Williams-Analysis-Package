function quantifyRemapping(folder)
%     dat = getTrialData(folder);
    load('combinedData');
%     allSim = [];
%     for mi = 1:length(dat)
%         for di = 1:length(dat(mi).maps.rooms(1,:))
%             doM = cat(4,dat(mi).maps.rooms{:,di});
%             sim = nan(3,3);
%             for i = 1:length(doM(1,1,1,:))
%                 for j = i+1:length(doM(1,1,1,:))
%                     m1 = doM(:,:,:,i);
%                     m2 = doM(:,:,:,j);
%                     isGood = ~isnan(m1)&~isnan(m2);
%                     sim(i,j) = corr(m1(isGood),m2(isGood));
%                 end
%             end
%             allSim = cat(3,allSim,sim);
%         end
%     end
%     
%     tmp = [];
%     for i = 1:length(allSim(1,1,:))
%         tmp2 = allSim(:,:,i);
%         tmp = [tmp; nanmean(tmp2(:))];
%     end
    
    for mi = 1:length(dat)
        figure(1)
        set(gcf,'position',[50 50 900 600])
        subplot(2,3,mi)
        h = cumHist(dat(mi).shr,[0:0.01:1]);
        hold on
        plot([0 1],[0 1],'linestyle','--','color','k','linewidth',1)
        xlabel('p-value')
        ylabel('Cumulative Proportion')
        axis equal
        axis square
        drawnow
        if mi == length(dat)
            legend(h,[{'Day 1'} {'Day 2'} {'Day 3'} {'Day 4'} {'Day 5'} {'Day 6'}]);
        end
    end

    allSim = [];
    for mi = 1:length(dat)
        for di = 1:length(dat(mi).maps.rooms(1,:))
            doM = cat(4,dat(mi).maps.rooms{:,di});
            sim = nan(3,3);
            for i = 1:length(doM(1,1,1,:))
                for j = i+1:length(doM(1,1,1,:))
                    m1 = doM(:,:,:,i);
                    m2 = doM(:,:,:,j);
                    isGood = ~isnan(m1)&~isnan(m2);
                    sim(i,j) = corr(m1(isGood),m2(isGood));
                end
            end
            allSim = cat(3,allSim,sim);
        end
    end
    
%     tmp = permute(nanmean(nanmean(allSim,1),2),[3 2 1]);
end