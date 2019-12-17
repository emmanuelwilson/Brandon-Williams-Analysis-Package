function [TopDist,AngleDiff,DistDiff] = TopologyClusterFinder(ms,angle,distance,maxdist,pass,mouse,context)
%% Will colour coat the spatial footprints based on preferred firing distance or angle.
% INPUT: 
%       - distanceMat: Matrix of preferred firing distances
%       - angleMat: Matrix of preferred firing angles
%       - ms.SFPs: Miniscope structure, must containt "SFPs" field
%       containing the SpatialFootPrints. 
% OUTPUT:
%       Produces colourcoated figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Émmanuel Wilson

% distance = 0;                                                          % Change to true if you wish to produce distance based figure
goodcells = ms.SFPs(:,:,pass);
% mouse = 'A2JH3';
% context = 'A';

%% Angle
pref = angle(pass);
dist = distance(pass);
sname = 'Angle';

for i = 1 : length(goodcells(1,1,:))
    mtemp = goodcells(:,:,i);
    maskThresh = prctile(mtemp(find(mtemp)),90);
    mind = find(mtemp>=maskThresh);    
    mtemp(mind) = 1;
    mtemp(find(mtemp~=1)) = 0;
    stats = regionprops(mtemp);
    centroids(i,:) = stats.Centroid;
end

if ~isempty(pass) && length(centroids(:,1)) > 1
    for i = 1 : length(centroids)-1
        for j = i+1 : length(centroids)
            TopDist((i-1)*(length(centroids))+j-1) = norm(centroids(i,:) - centroids(j,:));
            AngleDiff((i-1)*(length(centroids))+j-1) = abs(pref(i) - pref(j));
            DistDiff((i-1)*(length(centroids))+j-1) = abs(dist(i) - dist(j));
        end
    end
    AngleDiff(AngleDiff>180)  = 360 - AngleDiff(AngleDiff>180);
    
    figure(1)
    scatter(TopDist,AngleDiff)
    title([mouse 'Context' context ' Preferred Firing ' sname 'vs Topology Distance'],'FontSize',26)
    ylabel('Angle Difference (degrees)','FontSize',18)
    xlabel('Centroid Distance (pixels)','FontSize',18)
    ylim([0 180])
    
    saveas(gcf,['Pref', sname,'TopologyDifference', context,'_', mouse, '.fig'])
    saveas(gcf,['Pref', sname,'TopologyDifference', context, '_', mouse, '.eps'])
    
    %% Distance
    sname = 'Distance';
    
    figure(2)
    scatter(TopDist,DistDiff)
    title([mouse 'Context' context ' Preferred Firing ' sname 'vs Topology Distance'],'FontSize',26)
    ylabel('Distance Difference (cm)','FontSize',18)
    xlabel('Centroid Distance (pixels)','FontSize',18)
    ylim([0 maxdist])
    
    saveas(gcf,['Pref', sname,'TopologyDifference', context,'_', mouse, '.fig'])
    saveas(gcf,['Pref', sname,'TopologyDifference', context, '_', mouse, '.eps'])
else
    TopDist = [];
    AngleDiff = [];
    DistDiff = [];
end
end