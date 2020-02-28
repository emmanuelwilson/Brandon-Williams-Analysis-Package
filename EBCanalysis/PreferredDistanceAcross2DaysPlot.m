%% preferred firing distance plot generator
function [] = PreferredDistanceAcross2DaysPlot(DistanceAcross, session1,session2,mouse,save)

h1 = figure;
for i = 1 : length(DistanceAcross(:,1))
    plot(DistanceAcross(i,:))
    hold on
end
ylim([0 37.5])
xlim([0.5 2.5])
title([mouse, ': Context', session1, ' - ', 'Context', session2, ' Preferred Distance'],'FontSize',18)
xticks([1 2])
xticklabels({session1,session2})
ylabel('Distance (cm)','FontSize',14)

if save
    saveas(h1,['DistancePlot_',session1,'_',session2,'_', mouse,'.fig'])    
end

end