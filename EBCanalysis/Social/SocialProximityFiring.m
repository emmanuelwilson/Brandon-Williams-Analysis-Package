%% Will identify the middle of an object and monitor calcium activity within proximity of the object
function out = SocialProximityFiring(ms,tracking,frameMap,ebcstats, objectsizeRad,ROIrad, FieldSize)

maxX = max(ebcstats.QPW(:,1));
maxY = max(ebcstats.QPW(:,2));
minX = min(ebcstats.QPW(:,1));
minY = min(ebcstats.QPW(:,2));

maxXob = max(ebcstats.QPO(:,1));
maxYob = max(ebcstats.QPO(:,2));
minXob = min(ebcstats.QPO(:,1));
minYob = min(ebcstats.QPO(:,2));

pix2cm= (maxX-minX)/FieldSize;
objectPix = pix2cm*objectsizeRad;
ROIradPix = pix2cm*ROIrad;
tracking = tracking(frameMap,:);

objCenter = [maxXob-(maxXob-minXob)/2,maxYob-(maxYob-minYob)/2];
R = objectPix + ROIradPix;
L = linspace(0,2*pi);
xv = R*cos(L)';
yv = R*sin(L)';
totalOccupancy = sum(((tracking(:,1)-objCenter(1)).^2+(tracking(:,2)-objCenter(2)).^2<=R^2));

for i = 1 : length(ms.FiltTraces(1,:))
    tracktemp = tracking(ms.deconvolvedSig(:,i)>0,:);
    trackObjActive(i) = sum(((tracktemp(:,1)-objCenter(1)).^2+(tracktemp(:,2)-objCenter(2)).^2<=R^2));        
    trackActive(i) = sum(ms.deconvolvedSig(:,i)>0);
    tracktemp = [];
end

Objactivitypercent = trackObjActive./totalOccupancy;
totalactivitypercent = trackObjActive./trackActive;

out.objectCenter = objCenter;
out.Radius = R;
out.pix2cm = pix2cm;
out.TotalFramesObjOccupancy = totalOccupancy;
out.ActiveFramesObjOccupancy = trackObjActive';
out.ObjectOccPercentActiveProximity = Objactivitypercent';
out.PercentProximityActiveVsTotalActive = totalactivitypercent';

end