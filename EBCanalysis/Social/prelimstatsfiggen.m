%{
prelimSocialStats1.Sep_CupLeftFollow = cat(2,prelimSocialStats1.Sep_CupLeftFollow,prelimSocialStats2.Sep_CupLeftFollow);
prelimSocialStats1.Sep_CupRightFollow = cat(2,prelimSocialStats1.Sep_CupRightFollow,prelimSocialStats2.Sep_CupRightFollow);
prelimSocialStats1.Sep_MouseLeft2RightFollow = cat(2,prelimSocialStats1.Sep_MouseLeft2RightFollow,prelimSocialStats2.Sep_MouseLeft2RightFollow);
prelimSocialStats1.Sep_MouseRight2LeftFollow = cat(2,prelimSocialStats1.Sep_MouseRight2LeftFollow,prelimSocialStats2.Sep_MouseRight2LeftFollow);
prelimSocialStats1.Left1passnum = cat(2,prelimSocialStats1.Left1passnum,prelimSocialStats2.Left1passnum);
prelimSocialStats1.Left2passnum = cat(2,prelimSocialStats1.Left2passnum,prelimSocialStats2.Left2passnum);
prelimSocialStats1.Right1passnum = cat(2,prelimSocialStats1.Right1passnum,prelimSocialStats2.Right1passnum);
prelimSocialStats1.Right2passnum = cat(2,prelimSocialStats1.Right2passnum,prelimSocialStats2.Right2passnum);
%}

percentLeft2RightMouse = prelimSocialStats1.Sep_MouseLeft2RightFollow./prelimSocialStats1.Left1passnum;
percentRight2LefttMouse = prelimSocialStats1.Sep_MouseRight2LeftFollow./prelimSocialStats1.Right1passnum;
percentLeftCup = prelimSocialStats1.Sep_CupLeftFollow./prelimSocialStats1.Left1passnum;
percentRightCup = prelimSocialStats1.Sep_CupRightFollow./prelimSocialStats1.Right1passnum;

c = hsv(length(prelimSocialStats1.Left1passnum));
% colororder(c);
figure
hold on
for i = 1 : length(prelimSocialStats1.Left1passnum)    
    scatter(1,prelimSocialStats1.Left1passnum(i),[],c(i,:))
    scatter(2,prelimSocialStats1.Right1passnum(i),[],c(i,:))
    scatter(3,prelimSocialStats1.Left2passnum(i),[],c(i,:))
    scatter(4,prelimSocialStats1.Right2passnum(i),[],c(i,:))
end

xticks([1:4])
xticklabels({'Left1', 'Right1','Left2','Right2'})
xlim([0.5 4.5])
title('Number of Passed Cells')
ylabel('Number of Cells')

hold off

figure
hold on
for i = 1 : length(prelimSocialStats1.Left1passnum)
    scatter(1, percentLeft2RightMouse(i),[],c(i,:))
    scatter(2,percentRight2LefttMouse(i),[],c(i,:))
    scatter(3,percentLeftCup(i),[],c(i,:))
    scatter(4,percentRightCup(i),[],c(i,:))
end
xticks([1:4])
xticklabels({'Left to Right', 'Right to Left', 'Left to Left', 'Right to Right'})
xlim([0.5 4.5])
title('Proportion of Persistant Cells')
ylabel('Proportion of Cells (#Passed 1st&2nd half/# Passed in 1st half)')
