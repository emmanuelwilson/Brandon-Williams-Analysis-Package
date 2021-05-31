socialProximity = diagProximity;

leftind = 'diagProxLeft';
rightind = 'diagProxRight';

left = diagProxLeft;
right = diagProxRight;
lmins = diagminsLeft;
rmins = diagminsRight;


n = length(socialProximity.socialProxSessions(1,:));
m = 4;
s = [];
c = hsv(n);
for i = 1 : n
    stemp = [1:m ; i,i,i,i];
    s = cat(1,s,stemp);
end
exploringtimetotal = [];
exploringtime1 = [];
exploringtime2= [];
figure(1)
figure(2)
figure(3)
for i = 1 : n
    temp = nan(m,1);
    temp1 = nan(m,1);
    temp2 = nan(m,1);
    for j = 1 : m
        if ~isempty(socialProximity.socialProxSessions{j,i})
            temp(j) = socialProximity.socialProxSessions{j,i}.TotalFramesObjOccupancy;
            temp1(j) = socialProximity.socialProxSessions{j,i}.Object1FramesObjOccupancy;
            temp2(j) = socialProximity.socialProxSessions{j,i}.Object2FramesObjOccupancy;
        end
    end
    exploringtimetotal = cat(2,exploringtimetotal,temp);
    exploringtime1 = cat(2,exploringtime1,temp1);
    exploringtime2 = cat(2,exploringtime2,temp2);
    figure(1)
    plot(exploringtimetotal(:,i))    
    hold on
    figure(2)
    plot(exploringtime1(:,i))
    hold on
    figure(3)
    plot(exploringtime2(:,i))
    hold on
end

figure(1)
title('Total Time spent around Cups')
ylabel('Time (frames)')
xticks([1:4])
xticklabels({'Habituation','Trial1','Trial2','Trial3'})
ylim([0 max(max(exploringtimetotal))+200])
xlim([0.5 4.5])
colororder(c)
figure(2)
title('Total Time spent around Left Cup')
ylabel('Time (frames)')
xticks([1:4])
xticklabels({'Habituation','Trial1','Trial2','Trial3'})
ylim([0 max(max(exploringtimetotal))+200])
xlim([0.5 4.5])
colororder(c)
figure(3)
title('Total Time spent around Right Cup')
ylabel('Time (frames)')
xticks([1:4])
xticklabels({'Habituation','Trial1','Trial2','Trial3'})
ylim([0 max(max(exploringtimetotal))+200])
xlim([0.5 4.5])
colororder(c)

scatdist = linspace(0.8,1.2,n);

figure(4)
figure(5)
figure(6)
topval = 0;
topvaldiff = 0;
minvaldiff = 0;
toprate = 0;
topratediff = 0;
minratediff = 0;
for i = 1 : n
    for j = 1 : m
        temp = [];
        temp1 = [];
        temp2 = [];
        if ~isempty(socialProximity.socialProxSessions{j,i})
            temp = socialProximity.socialProxSessions{j,i}.ObjectOccPercentActiveProximity(find(socialProximity.exclude{j,i}.SFPs));
            temp1 = socialProximity.socialProxSessions{j,i}.FiringRateSignal_norm_Ob1(find(socialProximity.exclude{j,i}.SFPs))/socialProximity.socialProxSessions{j,i}.Object1FramesObjOccupancy;
            temp2 = socialProximity.socialProxSessions{j,i}.FiringRateSignal_norm_Ob2(find(socialProximity.exclude{j,i}.SFPs))/socialProximity.socialProxSessions{j,i}.Object2FramesObjOccupancy;
            
            figure(4)
            scatter(ones(length(temp),1)*(scatdist(i)+(j-1)),temp,[],c(i,:))
            hold on
            figure(5)
            scatter(ones(length(temp),1)*(scatdist(i)+(j-1)),temp1,[],c(i,:))
            hold on
            figure(6)
            scatter(ones(length(temp),1)*(scatdist(i)+(j-1)),temp2,[],c(i,:))
            hold on
                        
            tempFiringdiff = socialProximity.socialProxSessions{j,i}.FiringRateSignal_norm_Ob1(find(socialProximity.exclude{j,i}.SFPs))/socialProximity.socialProxSessions{j,i}.Object1FramesObjOccupancy - socialProximity.socialProxSessions{j,i}.FiringRateSignal_norm_Ob2(find(socialProximity.exclude{j,i}.SFPs))/socialProximity.socialProxSessions{j,i}.Object2FramesObjOccupancy;
            tempFiringdiff = -tempFiringdiff; %invert = becomes Object 2 - object 1
            tempmeandiff = mean(tempFiringdiff);
            
            figure(7)
            scatter(tempFiringdiff,ones(length(tempFiringdiff),1)*(scatdist(i)+(j-1)),[],c(i,:))
            hold on
            scatter(tempmeandiff,(scatdist(i)+(j-1)), 150, 'k', 'Marker', 'x')
            if exist(leftind) && exist(rightind) && j == 4
                if ~isempty(find(left == i))
                    figure(7)
                    ind = find(left == i);
                    scatter(min(tempFiringdiff) - 0.01, (scatdist(i)+(j-1)), 150, 'k', 'Marker', '*')
                    text(min(tempFiringdiff) - 0.02,(scatdist(i)+(j-1)), [num2str(lmins(ind)) 'mins'])
                    figure(5)
                    scatter((scatdist(i)+(j-1)),max(temp1) + 0.01, 150, 'k', 'Marker', '*')
                else
                    figure(7)
                    ind = find(right == i);
                    scatter( max(tempFiringdiff) + 0.01, (scatdist(i)+(j-1)) , 150, 'k', 'Marker', '*')
                    text(max(tempFiringdiff) + 0.015, (scatdist(i)+(j-1)) , [num2str(rmins(ind)) 'mins'])
                    figure(6)
                    scatter((scatdist(i)+(j-1)),max(temp2) + 0.01, 150, 'k', 'Marker', '*')
                end
            end
            
            frate_nocup = (SocialProximityAll.socialProxSessions{j,i}.FiringRateNorm_Sig(SocialProximityAll.exclude{j,i}.SFPs)-SocialProximityAll.socialProxSessions{j,i}.FiringRateSignal_norm_Ob2(SocialProximityAll.exclude{j,i}.SFPs)-SocialProximityAll.socialProxSessions{j,i}.FiringRateSignal_norm_Ob1(SocialProximityAll.exclude{j,i}.SFPs)) / (SocialProximityAll.socialProxSessions{j,i}.TotalFrames-SocialProximityAll.socialProxSessions{j,i}.TotalFramesObjOccupancy);
            frate_cup = (SocialProximityAll.socialProxSessions{j,i}.FiringRateSignal_norm_Ob2(SocialProximityAll.exclude{j,i}.SFPs)+SocialProximityAll.socialProxSessions{j,i}.FiringRateSignal_norm_Ob1(SocialProximityAll.exclude{j,i}.SFPs))/SocialProximityAll.socialProxSessions{j,i}.TotalFramesObjOccupancy;
            fratediff = frate_cup - frate_nocup;
            fratediffmean = mean(fratediff);
            figure(8)            
            subplot(2,2,j)
            plot([1, 2], [mean(frate_cup), mean(frate_nocup)],'Color',c(i,:))
            hold on
            
            figure(9)
            scatter(ones(length(fratediff),1)*(scatdist(i)+(j-1)),fratediff,[],c(i,:))
            hold on
            scatter((scatdist(i)+(j-1)), fratediffmean, 150, 'k', 'Marker', 'x')
            
            if max(temp) > topval
                topval = max(temp);
            end
            if max(tempFiringdiff) > topvaldiff
                topvaldiff = max(tempFiringdiff);
            end
            if min(tempFiringdiff) < minvaldiff
                minvaldiff = min(tempFiringdiff);
            end
            if (mean(frate_nocup)) > toprate || (mean(frate_cup)) > toprate
                if mean(frate_nocup) < mean(frate_cup)
                    toprate = mean(frate_cup);
                else
                    toprate = mean(frate_nocup);
                end
            end
            if max(fratediff) > topratediff
                topratediff = max(fratediff);
            end
            if min(fratediff) < minratediff
                minratediff = min(fratediff);
            end
        end
    end
end

figure(4)
title('Firing Rate around Cups')
ylabel('Object Activation (Active Frame/Proximity Frames)')
xticks([1:4])
xticklabels({'Habituation','Trial1','Trial2','Trial3'})
ylim([0 topval+0.05])
xlim([0.5 4.5])

figure(5)
title('Firing Rate around Left Cup')
ylabel('Object Activation (Active Frame/Proximity Frames)')
xticks([1:4])
xticklabels({'Habituation','Trial1','Trial2','Trial3'})
ylim([0 topval+0.05])
xlim([0.5 4.5])

figure(6)
title('Firing Rate around Right Cup')
ylabel('Object Activation (Active Frame/Proximity Frames)')
xticks([1:4])
xticklabels({'Habituation','Trial1','Trial2','Trial3'})
ylim([0 topval+0.05])
xlim([0.5 4.5])

figure(7)
title('Normalized Firing Rate Difference Between Cups')
xlabel('Difference in Firing Rate (Difference in Sum of normalized deconvolved Signal around Cups/Occupancy)')
yticks([1:4])
yticklabels({'Habituation','Trial1','Trial2','Trial3'})
xlim([minvaldiff-0.05 topvaldiff+0.05])
ylim([0.5 4.5])
vline(0, 'k')

figure(8)
ax1 = subplot(2,2,1);
title('Normalized Firing Rate Habituation')
xticks([1:2])
xticklabels({'Cups','No Cups'})
ylabel('(Normalized deconvolved Sum/Occupancy)')
xlim([0.5 2.5])
ylim([0 toprate+0.001])
ax2 = subplot(2,2,2);
title('Normalized Firing Rate Trial 1')
xticks([1:2])
xticklabels({'Cups','No Cups'})
ylabel('(Normalized deconvolved Sum/Occupancy)')
xlim([0.5 2.5])
ylim([0 toprate+0.001])
ax3 = subplot(2,2,3);
title('Normalized Firing Rate Trial 2')
xticks([1:2])
xticklabels({'Cups','No Cups'})
ylabel('(Normalized deconvolved Sum/Occupancy)')
xlim([0.5 2.5])
ylim([0 toprate+0.001])
ax4 = subplot(2,2,4);
title('Normalized Firing Rate Trial 3')
xticks([1:2])
xticklabels({'Cups','No Cups'})
ylabel('(Normalized deconvolved Sum/Occupancy)')
ylabel('Difference in Firing Rate (Sum of normalized deconvolved Signal/Occupancy)')
xticks([1:2])
xticklabels({'Cups','No Cups'})
ylim([0 toprate+0.001])
xlim([0.5 2.5])

figure(9)
title('Normalized Firing Rate Difference Cup vs non-Cup')
ylabel('Difference between sum of normalized Signal/Occupancy (cup - non-cup)')
xticks([1:4])
xticklabels({'Habituation','Trial1','Trial2','Trial3'})
ylim([minratediff-0.001 topratediff+0.001])
xlim([0.5 4.5])
hline(0, 'k')