%% Expirmental angle display between sessions
function [] = RingAngleRep(angles, session1, session2,mouse,save)
% angles = deg2rad(AngleDaysG23);
% session1 = 'A2';
% session2 = 'B';
% mouse = 'All Mice';

angles = deg2rad(angles + 90);

%generate ring
h1 = figure;
smallradius = 5; 
centerX = 0;
centerY = 0;
viscircles([centerX, centerY], smallradius,'Color','k');
axis square;
axis off

%overlay second ring
bigradius = 10; 
centerX = 0;
centerY = 0;
viscircles([centerX, centerY], bigradius,'Color','k');
axis square;
axis off
ylim([-bigradius-1,bigradius+1])
hline(0,'k')
vline(0,'k')

smallangles = angles(:,1);
bigangles = angles(:,2);
hold on
scatter(smallradius*cos(smallangles),smallradius*sin(smallangles))
scatter(bigradius*cos(bigangles),bigradius*sin(bigangles))

x = cat(2,smallradius*cos(smallangles),bigradius*cos(bigangles));
y = cat(2,smallradius*sin(smallangles),bigradius*sin(bigangles));
for i = 1 : length(angles(:,1))
    plot(x(i,:),y(i,:))
end

title([mouse,': Context ', session1, ' - Context ', session2,' Preferred Angle']) 
text(0.5,-1,['Context ',session1])
text(10.5*cos(11*pi/6),10.5*sin(11*pi/6),['Context ',session2])

text(0.5,bigradius+0.75,['0'])
text(0.5,-bigradius-0.75,['180'])
text(-bigradius-1,0.5,['90'])
text(bigradius+0.3,0.5,['270'])

if save
    saveas(h1,['AnlgeRing_',session1,'_',session2,'_', mouse,'.fig'])
    saveas(h1,['AnlgeRing_',session1,'_',session2,'_', mouse,'.eps'])
end
end