%% Will identify the middle of an object and monitor calcium activity within proximity of the object
function out = ObjectProximityFiring(ms,tracking,frameMap,Objectstats, objectsizeRad,ROIrad, FieldSizex,FieldSizey, badframes)
if ~isempty(badframes)
    frameMap(1:badframes) = [];
    ms.FiltTraces(1:badframes,:) = [];
    ms.deconvolvedSig(1:badframes,:) = [];
end
fps = 30;
mins = 10;
timeframes = mins*60;
lframes = timeframes*fps;
if length(frameMap) > lframes
    frameMap = frameMap(1:lframes);
    ms.FiltTraces = ms.FiltTraces(1:lframes,:);
    ms.deconvolvedSig = ms.deconvolvedSig(1:lframes,:);
end

maxX = max(Objectstats.QPW(:,1));
maxY = max(Objectstats.QPW(:,2));
minX = min(Objectstats.QPW(:,1));
minY = min(Objectstats.QPW(:,2));

pix2cmX= (maxX-minX)/FieldSizex;
pix2cmY= (maxY-minY)/FieldSizey;

objectPixX = pix2cmX*objectsizeRad;
objectPixY = pix2cmY*objectsizeRad;

ROIradPixX = pix2cmX*ROIrad;
ROIradPixY = pix2cmY*ROIrad;

tracking = tracking(frameMap,:);

objCenter1 = Objectstats.CentroidOb1;
objCenter2 = Objectstats.CentroidOb2;

R = objectPixX + ROIradPixX;
% Ry = objectPixY + ROIradPixY;

L = linspace(0,2*pi);
xv = R*cos(L)';
yv = R*sin(L)';

object1Occupancy = sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2)) - sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=objectPixX^2));
object2Occupancy = sum(((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2)) - sum(((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=objectPixX^2));
obj1OccFrames = ((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2) & ~((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=objectPixX^2);
obj2OccFrames  = ((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2) & ~((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=objectPixX^2);
objframtemp = cat(2,obj1OccFrames,obj2OccFrames);
objOccFrames = sum(objframtemp,2);
totalOccupancy = sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2)) + sum(((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2)) - sum(((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=objectPixX^2)) - sum(((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=objectPixX^2));
totalFrames = length(frameMap);

for i = 1 : length(ms.FiltTraces(1,:))    
    normCal = normalize(ms.FiltTraces(:,i),'range');
    normSig = normalize(ms.deconvolvedSig(:,i),'range');
    normSigtot(i,1) = sum(normSig);
    frO1_Cal(i,1) = sum(normCal((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2)) - sum(normCal((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=objectPixX^2));
    frO1_Sig(i,1) = sum(normSig((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=R^2)) - sum(normSig((tracking(:,1)-objCenter1(1)).^2+(tracking(:,2)-objCenter1(2)).^2<=objectPixX^2));
    frO2_Cal(i,1) = sum(normCal((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2)) - sum(normCal((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=objectPixX^2));
    frO2_Sig(i,1) = sum(normSig((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=R^2)) - sum(normSig((tracking(:,1)-objCenter2(1)).^2+(tracking(:,2)-objCenter2(2)).^2<=objectPixX^2));
    tracktemp = tracking(ms.deconvolvedSig(:,i)>0,:);
    trackObjActive(i) = sum(((tracktemp(:,1)-objCenter1(1)).^2+(tracktemp(:,2)-objCenter1(2)).^2<=R^2)) + sum(((tracktemp(:,1)-objCenter2(1)).^2+(tracktemp(:,2)-objCenter2(2)).^2<=R^2)) - sum(((tracktemp(:,1)-objCenter1(1)).^2+(tracktemp(:,2)-objCenter1(2)).^2<=objectPixX^2)) - sum(((tracktemp(:,1)-objCenter2(1)).^2+(tracktemp(:,2)-objCenter2(2)).^2<=objectPixX^2));
    trackActive(i) = sum(ms.deconvolvedSig(:,i)>0);
    tracktemp = [];
end

out.object1Center = objCenter1;
out.object2Center = objCenter2;
out.Radius = R;
out.pix2cm = pix2cmX;
out.Object1FramesObjOccupancy = object1Occupancy;
out.Object2FramesObjOccupancy = object2Occupancy;
out.TotalFramesObjOccupancy = totalOccupancy;
out.ActiveFramesObjOccupancy = trackObjActive';
out.FiringRateNorm_Sig = normSigtot;
out.FiringRateCalcium_norm_Ob1 = frO1_Cal;
out.FiringRateCalcium_norm_Ob2 = frO2_Cal;
out.FiringRateSignal_norm_Ob1 = frO1_Sig;
out.FiringRateSignal_norm_Ob2 = frO2_Sig;
out.TotalFrames = totalFrames;
out.ObjectRadius = objectsizeRad;
out.ProximityRadius = ROIrad;
out.FieldSize = [FieldSizex,FieldSizey];
out.Object1Occ = obj1OccFrames;
out.Object2Occ = obj2OccFrames;
out.ObjectOcc = objOccFrames;
end