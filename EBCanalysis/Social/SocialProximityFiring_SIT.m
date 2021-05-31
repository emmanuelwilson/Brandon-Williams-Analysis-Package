%% Will identify the middle of an object and monitor calcium activity within proximity of the object
function out = SocialProximityFiring_SIT(ms,tracking,frameMap,ebcstats, objectsizeRad,ROIrad, FieldSize, badframes)
if ~isempty(badframes)
    frameMap(1:badframes) = [];
    ms.FiltTraces(1:badframes,:) = [];
    ms.deconvolvedSig(1:badframes,:) = [];
end

maxX = max(ebcstats.QPW(:,1));
maxY = max(ebcstats.QPW(:,2));
minX = min(ebcstats.QPW(:,1));
minY = min(ebcstats.QPW(:,2));

maxXob1 = max(ebcstats.QPO(:,1));
maxYob1 = max(ebcstats.QPO(:,2));
minXob1 = min(ebcstats.QPO(:,1));
minYob1 = min(ebcstats.QPO(:,2));

pix2cm= (maxX-minX)/FieldSize;
objectPix = pix2cm*objectsizeRad;
ROIradPix = pix2cm*ROIrad;
tracking = tracking(frameMap,:);

objCenter1 = [maxXob1-(maxXob1-minXob1)/2,maxYob1-(maxYob1-minYob1)/2];
R = objectPix + ROIradPix;
L = linspace(0,2*pi);
xv = R*cos(L)';
yv = R*sin(L)';
object1Occupancy = sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2));
totalOccupancy = sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2));
totalFrames = length(frameMap);

for i = 1 : length(ms.FiltTraces(1,:))
    normCal = normalize(ms.FiltTraces(:,i),'range');
    normSig = normalize(ms.deconvolvedSig(:,i),'range');
    normSigtot(i,1) = sum(normSig);
    frO1_Cal(i,1) = sum(normCal((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2));
    frO1_Sig(i,1) = sum(normSig((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2));
    tracktemp = tracking(ms.deconvolvedSig(:,i)>0,:);
    trackObjActive(i) = sum(((tracktemp(:,1)-objCenter1(1)).^2+(tracktemp(:,2)-objCenter1(2)).^2<=R^2));
    trackActive(i) = sum(ms.deconvolvedSig(:,i)>0);
    trackObj1Active(i) = sum(((tracktemp(:,1)-objCenter1(1)).^2+(tracktemp(:,2)-objCenter1(2)).^2<=R^2));   
    tracktemp = [];
end

Obj1activitypercent = trackObj1Active./totalOccupancy; 
Objactivitypercent = trackObjActive./totalOccupancy;
totalactivitypercent = trackObjActive./trackActive;

out.object1Center = objCenter1;
out.Radius = R;
out.pix2cm = pix2cm;
out.QPW = ebcstats.QPW;
out.QPO = ebcstats.QPO;
out.Object1FramesObjOccupancy = object1Occupancy;
out.TotalFramesObjOccupancy = totalOccupancy;
out.ActiveFramesObjOccupancy = trackObjActive';
out.Object1OccPercentActiveProximity = Obj1activitypercent';
out.ObjectOccPercentActiveProximity = Objactivitypercent';
out.PercentProximityActiveVsTotalActive = totalactivitypercent';
out.FiringRateNorm_Sig = normSigtot;
out.FiringRateCalcium_norm_Ob1 = frO1_Cal;
out.FiringRateSignal_norm_Ob1 = frO1_Sig;
out.TotalFrames = totalFrames;
out.ObjectRadius = objectsizeRad;
out.ProximityRadius = ROIrad;
out.FieldSize = FieldSize;
end