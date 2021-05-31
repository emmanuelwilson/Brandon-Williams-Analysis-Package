objectProximity = ObjectProximityAll;

n = length(objectProximity.ObjectProxSessions(1,:));
m = 2;
s = [];
c = hsv(n);
for i = 1 : n
    stemp = [1:m ; i,i];
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
        if ~isempty(objectProximity.ObjectProxSessions{j,i})
            temp(j) = objectProximity.ObjectProxSessions{j,i}.TotalFramesObjOccupancy;
            temp1(j) = objectProximity.ObjectProxSessions{j,i}.Object1FramesObjOccupancy;
            temp2(j) = objectProximity.ObjectProxSessions{j,i}.Object2FramesObjOccupancy;
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
title('Total Time spent around Objects')
ylabel('Time (frames)')
xticks([1:2])
xticklabels({'Trial1','Trial2'})
ylim([0 max(max(exploringtimetotal))+200])
xlim([0.5 2.5])
colororder(c)
figure(2)
title('Total Time spent around Left Object')
ylabel('Time (frames)')
xticks([1:2])
xticklabels({'Trial1','Trial2'})
ylim([0 max(max(exploringtimetotal))+200])
xlim([0.5 2.5])
colororder(c)
figure(3)
title('Total Time spent around Right Object')
ylabel('Time (frames)')
xticks([1:2])
xticklabels({'Trial1','Trial2'})
ylim([0 max(max(exploringtimetotal))+200])
xlim([0.5 2.5])
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
        temp1 = [];
        temp2 = [];
        if ~isempty(objectProximity.ObjectProxSessions{j,i})            
            temp1 = objectProximity.ObjectProxSessions{j,i}.FiringRateSignal_norm_Ob1(find(objectProximity.exclude{j,i}.SFPs))/objectProximity.ObjectProxSessions{j,i}.Object1FramesObjOccupancy;
            temp2 = objectProximity.ObjectProxSessions{j,i}.FiringRateSignal_norm_Ob2(find(objectProximity.exclude{j,i}.SFPs))/objectProximity.ObjectProxSessions{j,i}.Object2FramesObjOccupancy;
                        
            figure(5)
            scatter(ones(length(temp1),1)*(scatdist(i)+(j-1)),temp1,[],c(i,:))
            hold on
            figure(6)
            scatter(ones(length(temp1),1)*(scatdist(i)+(j-1)),temp2,[],c(i,:))
            hold on
                        
            tempFiringdiff = objectProximity.ObjectProxSessions{j,i}.FiringRateSignal_norm_Ob1(find(objectProximity.exclude{j,i}.SFPs))/objectProximity.ObjectProxSessions{j,i}.Object1FramesObjOccupancy - objectProximity.ObjectProxSessions{j,i}.FiringRateSignal_norm_Ob2(find(objectProximity.exclude{j,i}.SFPs))/objectProximity.ObjectProxSessions{j,i}.Object2FramesObjOccupancy;
            tempFiringdiff = -tempFiringdiff; %invert = becomes Object 2 - object 1
            tempmeandiff = mean(tempFiringdiff);
            
            figure(7)
            scatter(tempFiringdiff,ones(length(tempFiringdiff),1)*(scatdist(i)+(j-1)),[],c(i,:))
            hold on
            scatter(tempmeandiff, (scatdist(i)+(j-1)),150, 'k', 'Marker', 'x')
            if exist('objProxLeft') && exist('objProxRight') && j ==2 && exist('objminsLeft') && exist('objminsRight')
                if ~isempty(find(objProxLeft == i))
                    figure(7)
                    scatter(min(tempFiringdiff) - 0.005, (scatdist(i)+(j-1)) , 150, 'k', 'Marker', '*')
                    text(min(tempFiringdiff) - 0.015,(scatdist(i)+(j-1)), [num2str(objminsLeft(find(objProxLeft == i))) 'mins'])
                    figure(5)
                    scatter((scatdist(i)+(j-1)), max(temp1) + 0.01 , 150, 'k', 'Marker', '*')
                else
                    figure(7)
                    scatter(max(tempFiringdiff) + 0.005, (scatdist(i)+(j-1)) , 150, 'k', 'Marker', '*')
                    text(max(tempFiringdiff) + 0.01, (scatdist(i)+(j-1)), [num2str(objminsRight(find(objProxRight == i))) 'mins'])
                    figure(6)
                    scatter((scatdist(i)+(j-1)), max(temp2) + 0.01 , 150, 'k', 'Marker', '*')
                end
            end
            
            frate_nocup = (ObjectProximityAll.ObjectProxSessions{j,i}.FiringRateNorm_Sig(ObjectProximityAll.exclude{j,i}.SFPs)-ObjectProximityAll.ObjectProxSessions{j,i}.FiringRateSignal_norm_Ob2(ObjectProximityAll.exclude{j,i}.SFPs)-ObjectProximityAll.ObjectProxSessions{j,i}.FiringRateSignal_norm_Ob1(ObjectProximityAll.exclude{j,i}.SFPs)) / (ObjectProximityAll.ObjectProxSessions{j,i}.TotalFrames-ObjectProximityAll.ObjectProxSessions{j,i}.TotalFramesObjOccupancy);
            frate_cup = (ObjectProximityAll.ObjectProxSessions{j,i}.FiringRateSignal_norm_Ob2(ObjectProximityAll.exclude{j,i}.SFPs)+ObjectProximityAll.ObjectProxSessions{j,i}.FiringRateSignal_norm_Ob1(ObjectProximityAll.exclude{j,i}.SFPs))/ObjectProximityAll.ObjectProxSessions{j,i}.TotalFramesObjOccupancy;
            fratediff = frate_cup - frate_nocup;
            fratediffmean = mean(fratediff);
            
            figure(4)
            scatter(ones(length(frate_cup),1)*(scatdist(i)+(j-1)),frate_cup,[],c(i,:))
            hold on
            
            figure(8)            
            subplot(1,2,j)
            plot([1, 2], [mean(frate_cup), mean(frate_nocup)],'Color',c(i,:))
            hold on
            
            figure(9)
            scatter(ones(length(fratediff),1)*(scatdist(i)+(j-1)),fratediff,[],c(i,:))
            hold on
            scatter((scatdist(i)+(j-1)), fratediffmean, 150, 'k', 'Marker', 'x')
            
            if max(frate_cup) > topval
                topval = max(frate_cup);
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
title('Firing Rate around Objects')
ylabel('Firing Rate (Normalized Deconvolved Signal/# of frames in Proximity)')
xticks([1:2])
xticklabels({'Trial1','Trial2'})
ylim([0 topval+0.05])
xlim([0.5 2.5])

figure(5)
title('Firing Rate around Left Object')
ylabel('Firing Rate (Normalized Deconvolved Signal/# of frames in Proximity)')
xticks([1:2])
xticklabels({'Trial1','Trial2'})
ylim([0 topval+0.05])
xlim([0.5 2.5])

figure(6)
title('Firing Rate around Right Object')
ylabel('Firing Rate (Normalized Deconvolved Signal/# of frames in Proximity)')
xticks([1:2])
xticklabels({'Trial1','Trial2'})
ylim([0 topval+0.05])
xlim([0.5 2.5])

figure(7)
title('Normalized Firing Rate Difference Between Objects')
xlabel('Difference in Firing Rate (Difference in Sum of normalized deconvolved Signal around Cups/Occupancy)')
yticks([1:2])
yticklabels({'Trial1','Trial2'})
xlim([minvaldiff-0.05 topvaldiff+0.05])
ylim([0.5 2.5])
vline(0, 'k')

figure(8)
ax1 = subplot(1,2,1);
title('Normalized Firing Rate Trial 1')
xticks([1:2])
xticklabels({'Objectss','No Objects'})
ylabel('(Normalized deconvolved Sum/Occupancy)')
xlim([0.5 2.5])
ylim([0 toprate+0.001])
ax2 = subplot(1,2,2);
title('Normalized Firing Rate Trial 2')
xticks([1:2])
xticklabels({'Objects','No Objects'})
ylabel('(Normalized deconvolved Sum/Occupancy)')
xlim([0.5 2.5])
ylim([0 toprate+0.001])

figure(9)
title('Normalized Firing Rate Difference Objects vs non-Objects')
ylabel('Difference between sum of normalized Signal/Occupancy (Objects- non-Objects)')
xticks([1:2])
xticklabels({'Trial1','Trial2'})
ylim([minratediff-0.001 topratediff+0.001])
xlim([0.5 2.5])
hline(0, 'k')