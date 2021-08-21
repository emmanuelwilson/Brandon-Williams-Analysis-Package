%% Will identify the middle of an object and monitor calcium activity within proximity of the object
function out = SocialProximityFiring_V3(ms,tracking,frameMap,ebcstats, objectsizeRad,ROIrad, FieldSize, badframes)
if ~isempty(badframes)
    frameMap(1:badframes) = [];
    ms.FiltTraces(1:badframes,:) = [];
    ms.deconvolvedSig(1:badframes,:) = [];
end
fps = 30;
timeframes = 5*60;
lframes = timeframes*fps;
if length(frameMap) > lframes
    frameMap = frameMap(1:lframes);
    ms.FiltTraces = ms.FiltTraces(1:lframes,:);
    ms.deconvolvedSig = ms.deconvolvedSig(1:lframes,:);
end

maxX = max(ebcstats.QPW(:,1));
maxY = max(ebcstats.QPW(:,2));
minX = min(ebcstats.QPW(:,1));
minY = min(ebcstats.QPW(:,2));

% maxXob1 = max(ebcstats.QPOL(:,1));
% maxYob1 = max(ebcstats.QPOL(:,2));
% minXob1 = min(ebcstats.QPOL(:,1));
% minYob1 = min(ebcstats.QPOL(:,2));
% 
% maxXob2 = max(ebcstats.QPOR(:,1));
% maxYob2 = max(ebcstats.QPOR(:,2));
% minXob2 = min(ebcstats.QPOR(:,1));
% minYob2 = min(ebcstats.QPOR(:,2));

pix2cm= (maxX-minX)/FieldSize;
objectPix = pix2cm*objectsizeRad;
ROIradPix = pix2cm*ROIrad;
tracking = tracking(frameMap,:);

objCenter1 = ebcstats.QPOL;% [maxXob1-(maxXob1-minXob1)/2,maxYob1-(maxYob1-minYob1)/2];
objCenter2 = ebcstats.QPOR;%[maxXob2-(maxXob2-minXob2)/2,maxYob2-(maxYob2-minYob2)/2];
R = objectPix + ROIradPix;
L = linspace(0,2*pi);
xv = R*cos(L)';
yv = R*sin(L)';
object1Occupancy = sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2));
object2Occupancy = sum(((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2));
totalOccupancy = sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2)) + sum(((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2));
totalFrames = length(frameMap);
obj1OccFrames = ((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2);
obj2OccFrames  = ((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2);
objframtemp = cat(2,obj1OccFrames,obj2OccFrames);
objOccFrames = sum(objframtemp,2);

for i = 1 : length(ms.FiltTraces(1,:))
    normCal = normalize(ms.FiltTraces(:,i),'range');
    normSig = normalize(ms.deconvolvedSig(:,i),'range');
    normSigtot(i,1) = sum(normSig);
    frO1_Cal(i,1) = sum(normCal((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2));
    frO1_Sig(i,1) = sum(normSig((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2));
    frO2_Cal(i,1) = sum(normCal((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2));
    frO2_Sig(i,1) = sum(normSig((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2));
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
out.Object1OccPercentActiveProximity = Obj1activitypercent';
out.Object2OccPercentActiveProximity = Obj2activitypercent';
out.ObjectOccPercentActiveProximity = Objactivitypercent';
out.PercentProximityActiveVsTotalActive = totalactivitypercent';
out.FiringRateNorm_Sig = normSigtot;
out.FiringRateCalcium_norm_Ob1 = frO1_Cal;
out.FiringRateCalcium_norm_Ob2 = frO2_Cal;
out.FiringRateSignal_norm_Ob1 = frO1_Sig;
out.FiringRateSignal_norm_Ob2 = frO2_Sig;
out.TotalFrames = totalFrames;
out.ObjectRadius = objectsizeRad;
out.ProximityRadius = ROIrad;
out.FieldSize = FieldSize;
out.Object1Occ = obj1OccFrames;
out.Object2Occ = obj2OccFrames;
out.ObjectOcc = objOccFrames;
end