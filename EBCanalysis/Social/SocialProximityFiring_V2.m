%% Will identify the middle of an object and monitor calcium activity within proximity of the object
function out = SocialProximityFiring_V2(ms,tracking,frameMap,ebcstats_Ob1,ebcstats_Ob2, objectsizeRad,ROIrad, FieldSize)

maxX = max(ebcstats_Ob1.QPW(:,1));
maxY = max(ebcstats_Ob1.QPW(:,2));
minX = min(ebcstats_Ob1.QPW(:,1));
minY = min(ebcstats_Ob1.QPW(:,2));

maxXob1 = max(ebcstats_Ob1.QPO(:,1));
maxYob1 = max(ebcstats_Ob1.QPO(:,2));
minXob1 = min(ebcstats_Ob1.QPO(:,1));
minYob1 = min(ebcstats_Ob1.QPO(:,2));

maxXob2 = max(ebcstats_Ob2.QPO(:,1));
maxYob2 = max(ebcstats_Ob2.QPO(:,2));
minXob2 = min(ebcstats_Ob2.QPO(:,1));
minYob2 = min(ebcstats_Ob2.QPO(:,2));

pix2cm= (maxX-minX)/FieldSize;
objectPix = pix2cm*objectsizeRad;
ROIradPix = pix2cm*ROIrad;
tracking = tracking(frameMap,:);

objCenter1 = [maxXob1-(maxXob1-minXob1)/2,maxYob1-(maxYob1-minYob1)/2];
objCenter2 = [maxXob2-(maxXob2-minXob2)/2,maxYob2-(maxYob2-minYob2)/2];
R = objectPix + ROIradPix;
L = linspace(0,2*pi);
xv = R*cos(L)';
yv = R*sin(L)';
object1Occupancy = sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2));
object2Occupancy = sum(((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2));
totalOccupancy = sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2)) + sum(((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2));

for i = 1 : length(ms.FiltTraces(1,:))
    tracktemp = tracking(ms.deconvolvedSig(:,i)>0,:);
    trackObjActive(i) = sum(((tracktemp(:,1)-objCenter1(1)).^2+(tracktemp(:,2)-objCenter1(2)).^2<=R^2)) + sum(((tracktemp(:,1)-objCenter2(1)).^2+(tracktemp(:,2)-objCenter2(2)).^2<=R^2));
    trackActive(i) = sum(ms.deconvolvedSig(:,i)>0);
    trackObj1Active(i) = sum(((tracktemp(:,1)-objCenter1(1)).^2+(tracktemp(:,2)-objCenter1(2)).^2<=R^2));
    trackObj2Active(i) = sum(((tracktemp(:,1)-objCenter2(1)).^2+(tracktemp(:,2)-objCenter2(2)).^2<=R^2));
    tracktemp = [];
end

Obj1activitypercent = trackObj1Active./totalOccupancy; 
Obj2activitypercent = trackObj2Active./totalOccupancy;
Objactivitypercent = trackObjActive./totalOccupancy;
totalactivitypercent = trackObjActive./trackActive;

out.object1Center = objCenter1;
out.object2Center = objCenter2;
out.Radius = R;
out.pix2cm = pix2cm;
out.Object1FramesObjOccupancy = object1Occupancy;
out.Object2FramesObjOccupancy = object2Occupancy;
out.TotalFramesObjOccupancy = totalOccupancy;
out.ActiveFramesObjOccupancy = trackObjActive';
out.Object1OccPercentActiveProximity = Obj1activitypercent;
out.Object2OccPercentActiveProximity = Obj2activitypercent;
out.ObjectOccPercentActiveProximity = Objactivitypercent';
out.PercentProximityActiveVsTotalActive = totalactivitypercent';
out.ObjectRadius = objectsizeRad;
out.ProximityRadius = ROIrad;
out.FieldSize = FieldSize;
end