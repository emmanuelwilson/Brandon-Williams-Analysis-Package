scatmatSubtot = scat8;
endstep = 0.013;
stepSize = 0.0005;
fstep = 0.0013;
xbin = [fstep:stepSize:endstep];

vals = nan(length(xbin)-1,2);
for i = 1:length(xbin)-1
%     test = scatmattot(scatmattB(:,1)>xbin(i) & scatmattotB(:,1) < xbin(i+1));
    vals(i,1) = nanmean(scatmatSubtot(scatmatSubtot(:,2)>xbin(i) & scatmatSubtot(:,2) < xbin(i+1),1));    
    if xbin(i)> 0.0838
        length(find(scatmatSubtot(scatmatSubtot(:,2)>xbin(i) & scatmatSubtot(:,2) < xbin(i+1),1)));
    end
    vals(i,2) = nanstd(scatmatSubtot(scatmatSubtot(:,2)>xbin(i) & scatmatSubtot(:,2) < xbin(i+1),1));
end

maxMF = 0.015;  %Max Mean Firing (active frames/total frames)  %0.0905
minMF = 0.0;
stepSize = 0.00001;
xbin2 = [fstep:stepSize:endstep];
longvals = spline(xbin(1:end-1),vals(:,1),xbin2);
%Get rid of values that are interpolated past the max value (will vary)
longvals(xbin2>maxMF | xbin2<minMF) = [];
xbin2(xbin2>maxMF | xbin2<minMF) = [];
%plot interpolated data
figure
plot(xbin2,longvals)
ylabel('Mean Shuffled MRL')
xlabel('Mean Firing Rate')
%find the 2nd order changepoint
findchangepts(longvals,'MaxNumChanges',1)
ylabel('Mean Shuffled MRL')
xlabel('Time Bins')

figure
hold on
plot(xbin(1:end-1)+stepSize./2,vals(:,1), 'color', 'r');
plot(xbin(1:end-1)+stepSize./2,vals(:,1)+3.*vals(:,2), 'color', 'k');
plot(xbin(1:end-1)+stepSize./2,vals(:,1)-3.*vals(:,2), 'color', 'k');
ylabel('Mean Shuffled MRL')
xlabel('Mean Firing Rate')