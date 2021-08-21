%% Will identify the middle of an object and monitor calcium activity within proximity of the object
function out = SocialProximityFiring_SIT_V2(ms,tracking,frameMap,ebcstats,objstatshab, objectsizeRad,ROIrad, FieldSize, badframes, startframe)

frameMaptrial = frameMap;
frameMaphab = frameMap;
mshab = ms;
frameMaptrial(1:startframe) = [];
% ms.FiltTraces(1:startframe,:) = [];
% ms.deconvolvedSig(1:startframe,:) = [];

frameMaphab(1:badframes) = [];
% ms.FiltTraces(1:badframes,:) = [];
% ms.deconvolvedSig(1:badframes,:) = [];

fps = 30;
timeframes = 5*60;
lframes = timeframes*fps;

if lframes + badframes < startframe
    frameMaphab = frameMaphab(1:lframes);
    mshab.FiltTraces = ms.FiltTraces(badframes+1:badframes+lframes,:);
    mshab.deconvolvedSig = ms.deconvolvedSig(badframes+1:badframes+lframes,:);
else
    frameMaphab = frameMaphab(1:startframe);
    mshab.FiltTraces = ms.FiltTraces(badframes+1:badframes+startframe,:);
    mshab.deconvolvedSig = ms.deconvolvedSig(badframes+1:badframes+startframe,:);
end

if length(frameMaptrial) > lframes
    frameMaptrial = frameMaptrial(1:lframes);
    ms.FiltTraces = ms.FiltTraces(startframe+1:startframe+lframes,:);
    ms.deconvolvedSig = ms.deconvolvedSig(startframe+1:startframe+lframes,:);
else    
    ms.FiltTraces = ms.FiltTraces(startframe+1:startframe+length(frameMaptrial),:);
    ms.deconvolvedSig = ms.deconvolvedSig(startframe+1:startframe+length(frameMaptrial),:);
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
trackingtrial = tracking(frameMaptrial,:);
trackinghab = tracking(frameMaphab,:);

objCenter1 = [maxXob1-(maxXob1-minXob1)/2,maxYob1-(maxYob1-minYob1)/2];
R = objectPix + ROIradPix;
L = linspace(0,2*pi);
xv = R*cos(L)';
yv = R*sin(L)';
object1Occupancy = sum(((trackingtrial(:,1)-objCenter1(1)).^2+(trackingtrial(:,2)-objCenter1(2)).^2<=R^2));
totalOccupancy = sum(((trackingtrial(:,1)-objCenter1(1)).^2+(trackingtrial(:,2)-objCenter1(2)).^2<=R^2));
totalFrames = length(frameMaptrial);
obj1OccFrames = ((trackingtrial(:,1)-objCenter1(1)).^2+(trackingtrial(:,2)-objCenter1(2)).^2<=R^2);

object1Occupancyhab = sum(((trackinghab(:,1)-objstatshab.CentroidOb1(1)).^2+(trackinghab(:,2)-objstatshab.CentroidOb1(2)).^2<=R^2));
totalFrameshab = length(frameMaphab);
obj1OccFramesHab = ((trackinghab(:,1)-objstatshab.CentroidOb1(1)).^2+(trackinghab(:,2)-objstatshab.CentroidOb1(2)).^2<=R^2);


for i = 1 : length(ms.FiltTraces(1,:))
    normCal = normalize(ms.FiltTraces(:,i),'range');
    normSig = normalize(ms.deconvolvedSig(:,i),'range');
    normSigtot(i,1) = sum(normSig);
    frO1_Cal(i,1) = sum(normCal((trackingtrial(:,1)-objCenter1(1)).^2+(trackingtrial(:,2)-objCenter1(2)).^2<=R^2));
    frO1_Sig(i,1) = sum(normSig((trackingtrial(:,1)-objCenter1(1)).^2+(trackingtrial(:,2)-objCenter1(2)).^2<=R^2));
    tracktemp = trackingtrial(ms.deconvolvedSig(:,i)>0,:);
    trackObjActive(i) = sum(((tracktemp(:,1)-objCenter1(1)).^2+(tracktemp(:,2)-objCenter1(2)).^2<=R^2));
    trackActive(i) = sum(ms.deconvolvedSig(:,i)>0);
    trackObj1Active(i) = sum(((tracktemp(:,1)-objCenter1(1)).^2+(tracktemp(:,2)-objCenter1(2)).^2<=R^2));   
    tracktemp = [];
    
    normCalHab = normalize(mshab.FiltTraces(:,i),'range');
    normSigHab = normalize(mshab.deconvolvedSig(:,i),'range');
    normSigtotHab(i,1) = sum(normSigHab);
    frO1_CalHab(i,1) = sum(normCalHab((trackinghab(:,1)-objstatshab.CentroidOb1(1)).^2+(trackinghab(:,2)-objstatshab.CentroidOb1(2)).^2<=R^2));
    frO1_SigHab(i,1) = sum(normSigHab((trackinghab(:,1)-objstatshab.CentroidOb1(1)).^2+(trackinghab(:,2)-objstatshab.CentroidOb1(2)).^2<=R^2));
    tracktempHab = trackinghab(mshab.deconvolvedSig(:,i)>0,:);
    trackObjActiveHab(i) = sum(((tracktempHab(:,1)-objstatshab.CentroidOb1(1)).^2+(tracktempHab(:,2)-objstatshab.CentroidOb1(2)).^2<=R^2));
    trackActiveHab(i) = sum(mshab.deconvolvedSig(:,i)>0);
    trackObj1ActiveHab(i) = sum(((tracktempHab(:,1)-objstatshab.CentroidOb1(1)).^2+(tracktempHab(:,2)-objstatshab.CentroidOb1(2)).^2<=R^2));   
    tracktempHab = [];
end

Obj1activitypercent = trackObj1Active./totalOccupancy; 
Objactivitypercent = trackObjActive./totalOccupancy;
totalactivitypercent = trackObjActive./trackActive;

out.object1Center = objCenter1;
out.objectCenterHab = objstatshab.CentroidOb1;
out.Radius = R;
out.pix2cm = pix2cm;
out.QPW = ebcstats.QPW;
out.QPO = ebcstats.QPO;
out.Object1FramesObjOccupancy = object1Occupancy;
out.ObjectHabFramesObjOccupancy = object1Occupancyhab;
out.TotalFramesObjOccupancy = totalOccupancy;
out.ActiveFramesObjOccupancy = trackObjActive';
out.Object1OccPercentActiveProximity = Obj1activitypercent';
out.ObjectOccPercentActiveProximity = Objactivitypercent';
out.PercentProximityActiveVsTotalActive = totalactivitypercent';
out.FiringRateNorm_Sig = normSigtot;
out.FiringRateCalcium_norm_Ob1 = frO1_Cal;
out.FiringRateSignal_norm_Ob1 = frO1_Sig;
out.FiringRateCalcium_norm_Ob1Hab = frO1_CalHab;
out.FiringRateSignal_norm_Ob1Hab = frO1_SigHab;
out.TotalFrames = totalFrames;
out.TotalFramesHab = totalFrameshab;
out.ObjectRadius = objectsizeRad;
out.ProximityRadius = ROIrad;
out.FieldSize = FieldSize;
out.Object1Occ = obj1OccFrames;
out.ObjectOccHab = obj1OccFramesHab;
end