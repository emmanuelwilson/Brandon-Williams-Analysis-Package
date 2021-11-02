%% Create figure for concactenated Sessions SIT data


secframes2 = 60;
secframes3 = 90;

load('SITproximityData.mat')
load('SITnovelty_mins.mat')
cdold = pwd;
mousecount = 1;

interactionlength1 = [];
interactionlength2 = [];
interactionlengthHab1 = [];
interactionlengthHab2 = [];

chronoall = cell(length(SITnovelty_mins(:,1)),4);

ChronoAllcells120Hab1Same = [];
ChronoAllcells120Trial1Same = [];
ChronoAllcells120Hab2Same = [];
ChronoAllcells120Trial2Same = [];
ChronoAllcells120Hab1Diff = [];
ChronoAllcells120Trial1Diff = [];
ChronoAllcells120Hab2Diff = [];
ChronoAllcells120Trial2Diff = [];
ChronoAllcells30Hab1Same = [];
ChronoAllcells30Trial1Same = [];
ChronoAllcells30Hab2Same = [];
ChronoAllcells30Trial2Same = [];
ChronoAllcells30Hab1Diff = [];
ChronoAllcells30Trial1Diff = [];
ChronoAllcells30Hab2Diff = [];
ChronoAllcells30Trial2Diff = [];

MeanFiringAllcells120Trial1Same = [];
MeanFiringAllcells120Trial2Same = [];
MeanFiringAllcells120Trial1Diff = [];
MeanFiringAllcells120Trial2Diff = [];
MeanFiringAllcells30Trial1Same = [];
MeanFiringAllcells30Trial2Same = [];
MeanFiringAllcells30Trial1Diff = [];
MeanFiringAllcells30Trial2Diff = [];
MeanFiringAllcells120Trial1SameRaw = [];
MeanFiringAllcells120Trial2SameRaw = [];
MeanFiringAllcells120Trial1DiffRaw = [];
MeanFiringAllcells120Trial2DiffRaw = [];
MeanFiringAllcells30Trial1SameRaw = [];
MeanFiringAllcells30Trial2SameRaw = [];
MeanFiringAllcells30Trial1DiffRaw = [];
MeanFiringAllcells30Trial2DiffRaw = [];
MeanFiringAllcells120Trial1Same_passed = [];
MeanFiringAllcells120Trial2Same_passed = [];
MeanFiringAllcells120Trial1Diff_passed = [];
MeanFiringAllcells120Trial2Diff_passed = [];
MeanFiringAllcells30Trial1Same_passed = [];
MeanFiringAllcells30Trial2Same_passed = [];
MeanFiringAllcells30Trial1Diff_passed = [];
MeanFiringAllcells30Trial2Diff_passed = [];

MeanFiringAllcells120Hab1Same = [];
MeanFiringAllcells120Hab2Same = [];
MeanFiringAllcells120Hab1Diff = [];
MeanFiringAllcells120Hab2Diff = [];
MeanFiringAllcells30Hab1Same = [];
MeanFiringAllcells30Hab2Same = [];
MeanFiringAllcells30Hab1Diff = [];
MeanFiringAllcells30Hab2Diff = [];
MeanFiringAllcells120Hab1Same_passed = [];
MeanFiringAllcells120Hab2Same_passed = [];
MeanFiringAllcells120Hab1Diff_passed = [];
MeanFiringAllcells120Hab2Diff_passed = [];
MeanFiringAllcells30Hab1Same_passed = [];
MeanFiringAllcells30Hab2Same_passed = [];
MeanFiringAllcells30Hab1Diff_passed = [];
MeanFiringAllcells30Hab2Diff_passed = [];
MeanFiringAllcells120Hab1SameRaw = [];
MeanFiringAllcells120Hab2SameRaw = [];
MeanFiringAllcells120Hab1DiffRaw = [];
MeanFiringAllcells120Hab2DiffRaw = [];
MeanFiringAllcells30Hab1SameRaw = [];
MeanFiringAllcells30Hab2SameRaw = [];
MeanFiringAllcells30Hab1DiffRaw = [];
MeanFiringAllcells30Hab2DiffRaw = [];

MeanFiringAllcells120Trial1Same_passed_CR = [];
MeanFiringAllcells120Trial2Same_passed_CR = [];
MeanFiringAllcells120Trial1Diff_passed_CR = [];
MeanFiringAllcells120Trial2Diff_passed_CR = [];
MeanFiringAllcells30Trial1Same_passed_CR = [];
MeanFiringAllcells30Trial2Same_passed_CR = [];
MeanFiringAllcells30Trial1Diff_passed_CR = [];
MeanFiringAllcells30Trial2Diff_passed_CR = [];
MeanFiringAllcells120Hab1Same_passed_CR = [];
MeanFiringAllcells120Hab2Same_passed_CR = [];
MeanFiringAllcells120Hab1Diff_passed_CR = [];
MeanFiringAllcells120Hab2Diff_passed_CR = [];
MeanFiringAllcells30Hab1Same_passed_CR = [];
MeanFiringAllcells30Hab2Same_passed_CR = [];
MeanFiringAllcells30Hab1Diff_passed_CR = [];
MeanFiringAllcells30Hab2Diff_passed_CR = [];

%Trial1 Same vs Trial1 Diff
MeanFiringCellReg_30_S1D1 = [];
MeanFiringCellReg_30_Sn1D1 = [];
MeanFiringCellReg_30_S1Dn1 = [];
MeanFiringCellReg_30_Sn1Dn1 = [];
%Trial2 Same vs Trial1 Diff
MeanFiringCellReg_30_S2D1 = [];
MeanFiringCellReg_30_Sn2D1 = [];
MeanFiringCellReg_30_S2Dn1 = [];
MeanFiringCellReg_30_Sn2Dn1 = [];
%Trial2 Same vs Trial2 Diff
MeanFiringCellReg_30_S2D2 = [];
MeanFiringCellReg_30_Sn2D2 = [];
MeanFiringCellReg_30_S2Dn2 = [];
MeanFiringCellReg_30_Sn2Dn2 = [];
%Trial1 Same vs Trial2 Diff
MeanFiringCellReg_30_S1D2 = [];
MeanFiringCellReg_30_Sn1D2 = [];
MeanFiringCellReg_30_S1Dn2 = [];
MeanFiringCellReg_30_Sn1Dn2 = [];
%Hab1 Same vs Hab1 Diff
MeanFiringCellReg_30_Sh1Dh1 = [];
MeanFiringCellReg_30_Snh1Dh1 = [];
MeanFiringCellReg_30_Sh1Dnh1 = [];
MeanFiringCellReg_30_Snh1Dnh1 = [];
%Hab1 Same vs Hab2 Diff
MeanFiringCellReg_30_Sh1Dh2 = [];
MeanFiringCellReg_30_Snh1Dh2 = [];
MeanFiringCellReg_30_Sh1Dnh2 = [];
MeanFiringCellReg_30_Snh1Dnh2 = [];
%Hab2 Same vs Hab2 Diff
MeanFiringCellReg_30_Sh2Dh2 = [];
MeanFiringCellReg_30_Snh2Dh2 = [];
MeanFiringCellReg_30_Sh2Dnh2 = [];
MeanFiringCellReg_30_Snh2Dnh2 = [];
%Hab2 Same vs Hab1 Diff
MeanFiringCellReg_30_Sh2Dh1 = [];
MeanFiringCellReg_30_Snh2Dh1 = [];
MeanFiringCellReg_30_Sh2Dnh1 = [];
MeanFiringCellReg_30_Snh2Dnh1 = [];
%Trial1 Same vs Trial2 Same
MeanFiringCellReg_30_S1S2 = [];
MeanFiringCellReg_30_Sn1S2 = [];
MeanFiringCellReg_30_S1Sn2 = [];
MeanFiringCellReg_30_Sn1Sn2 = [];
%Trial1 Diff vs Trial2 Diff
MeanFiringCellReg_30_D1D2 = [];
MeanFiringCellReg_30_Dn1D2 = [];
MeanFiringCellReg_30_D1Dn2 = [];
MeanFiringCellReg_30_Dn1Dn2 = [];
%Trial1 Same vs Trial1 Diff
MeanFiringCellReg_120_S1D1 = [];
MeanFiringCellReg_120_Sn1D1 = [];
MeanFiringCellReg_120_S1Dn1 = [];
MeanFiringCellReg_120_Sn1Dn1 = [];
%Trial2 Same vs Trial1 Diff
MeanFiringCellReg_120_S2D1 = [];
MeanFiringCellReg_120_Sn2D1 = [];
MeanFiringCellReg_120_S2Dn1 = [];
MeanFiringCellReg_120_Sn2Dn1 = [];
%Trial2 Same vs Trial2 Diff
MeanFiringCellReg_120_S2D2 = [];
MeanFiringCellReg_120_Sn2D2 = [];
MeanFiringCellReg_120_S2Dn2 = [];
MeanFiringCellReg_120_Sn2Dn2 = [];
%Trial1 Same vs Trial2 Diff
MeanFiringCellReg_120_S1D2 = [];
MeanFiringCellReg_120_Sn1D2 = [];
MeanFiringCellReg_120_S1Dn2 = [];
MeanFiringCellReg_120_Sn1Dn2 = [];
%Hab1 Same vs Hab1 Diff
MeanFiringCellReg_120_Sh1Dh1 = [];
MeanFiringCellReg_120_Snh1Dh1 = [];
MeanFiringCellReg_120_Sh1Dnh1 = [];
MeanFiringCellReg_120_Snh1Dnh1 = [];
%Hab1 Same vs Hab2 Diff
MeanFiringCellReg_120_Sh1Dh2 = [];
MeanFiringCellReg_120_Snh1Dh2 = [];
MeanFiringCellReg_120_Sh1Dnh2 = [];
MeanFiringCellReg_120_Snh1Dnh2 = [];
%Hab2 Same vs Hab2 Diff
MeanFiringCellReg_120_Sh2Dh2 = [];
MeanFiringCellReg_120_Snh2Dh2 = [];
MeanFiringCellReg_120_Sh2Dnh2 = [];
MeanFiringCellReg_120_Snh2Dnh2 = [];
%Hab2 Same vs Hab1 Diff
MeanFiringCellReg_120_Sh2Dh1 = [];
MeanFiringCellReg_120_Snh2Dh1 = [];
MeanFiringCellReg_120_Sh2Dnh1 = [];
MeanFiringCellReg_120_Snh2Dnh1 = [];
%Trial1 Same vs Trial2 Same
MeanFiringCellReg_120_S1S2 = [];
MeanFiringCellReg_120_Sn1S2 = [];
MeanFiringCellReg_120_S1Sn2 = [];
MeanFiringCellReg_120_Sn1Sn2 = [];
%Trial1 Diff vs Trial2 Diff
MeanFiringCellReg_120_D1D2 = [];
MeanFiringCellReg_120_Dn1D2 = [];
MeanFiringCellReg_120_D1Dn2 = [];
MeanFiringCellReg_120_Dn1Dn2 = [];

%Colour Coated Firing Rate non-CellReg
%Trial1 Same vs Trial2 Same
MeanFiring_120_S1S2 = [];
MeanFiring_120_Sn1S2 = [];
MeanFiring_120_S1Sn2 = [];
MeanFiring_120_Sn1Sn2 = [];
%Trial1 Diff vs Trial2 Diff
MeanFiring_120_D1D2 = [];
MeanFiring_120_Dn1D2 = [];
MeanFiring_120_D1Dn2 = [];
MeanFiring_120_Dn1Dn2 = [];
%Hab1 Same vs Trial2 Same
MeanFiring_120_Sh1S2 = [];
MeanFiring_120_Snh1S2 = [];
MeanFiring_120_Sh1Sn2 = [];
MeanFiring_120_Snh1Sn2 = [];
%Hab1 Diff vs Trial2 Diff
MeanFiring_120_Dh1D2 = [];
MeanFiring_120_Dnh1D2 = [];
MeanFiring_120_Dh1Dn2 = [];
MeanFiring_120_Dnh1Dn2 = [];
%Trial1 Same vs Hab2 Same
MeanFiring_120_S1Sh2 = [];
MeanFiring_120_Sn1Sh2 = [];
MeanFiring_120_S1Snh2 = [];
MeanFiring_120_Sn1Snh2 = [];
%Trial1 Diff vs Hab2 Diff
MeanFiring_120_D1Dh2 = [];
MeanFiring_120_Dn1Dh2 = [];
MeanFiring_120_D1Dnh2 = [];
MeanFiring_120_Dn1Dnh2 = [];
%Hab1 Same vs Hab2 Same
MeanFiring_120_Sh1Sh2 = [];
MeanFiring_120_Snh1Sh2 = [];
MeanFiring_120_Sh1Snh2 = [];
MeanFiring_120_Snh1Snh2 = [];
%Hab1 Diff vs Hab2 Diff
MeanFiring_120_Dh1Dh2 = [];
MeanFiring_120_Dnh1Dh2 = [];
MeanFiring_120_Dh1Dnh2 = [];
MeanFiring_120_Dnh1Dnh2 = [];

%Trial1 Same vs Trial2 Same
MeanFiring_30_S1S2 = [];
MeanFiring_30_Sn1S2 = [];
MeanFiring_30_S1Sn2 = [];
MeanFiring_30_Sn1Sn2 = [];
%Trial1 Diff vs Trial2 Diff
MeanFiring_30_D1D2 = [];
MeanFiring_30_Dn1D2 = [];
MeanFiring_30_D1Dn2 = [];
MeanFiring_30_Dn1Dn2 = [];
%Hab1 Same vs Trial2 Same
MeanFiring_30_Sh1S2 = [];
MeanFiring_30_Snh1S2 = [];
MeanFiring_30_Sh1Sn2 = [];
MeanFiring_30_Snh1Sn2 = [];
%Hab1 Diff vs Trial2 Diff
MeanFiring_30_Dh1D2 = [];
MeanFiring_30_Dnh1D2 = [];
MeanFiring_30_Dh1Dn2 = [];
MeanFiring_30_Dnh1Dn2 = [];
%Trial1 Same vs Hab2 Same
MeanFiring_30_S1Sh2 = [];
MeanFiring_30_Sn1Sh2 = [];
MeanFiring_30_S1Snh2 = [];
MeanFiring_30_Sn1Snh2 = [];
%Trial1 Diff vs Hab2 Diff
MeanFiring_30_D1Dh2 = [];
MeanFiring_30_Dn1Dh2 = [];
MeanFiring_30_D1Dnh2 = [];
MeanFiring_30_Dn1Dnh2 = [];
%Hab1 Same vs Hab2 Same
MeanFiring_30_Sh1Sh2 = [];
MeanFiring_30_Snh1Sh2 = [];
MeanFiring_30_Sh1Snh2 = [];
MeanFiring_30_Snh1Snh2 = [];
%Hab1 Diff vs Hab2 Diff
MeanFiring_30_Dh1Dh2 = [];
MeanFiring_30_Dnh1Dh2 = [];
MeanFiring_30_Dh1Dnh2 = [];
MeanFiring_30_Dnh1Dnh2 = [];


Chrono120Hab1Same = [];
Chrono1201Same = [];
Chrono120Hab2Same = [];
Chrono1202Same = [];
Chrono120Hab1Diff = [];
Chrono1201Diff = [];
Chrono120Hab2Diff = [];
Chrono1202Diff = [];
Chrono30Hab1Same = [];
Chrono301Same = [];
Chrono30Hab2Same = [];
Chrono302Same = [];
Chrono30Hab1Diff = [];
Chrono301Diff = [];
Chrono30Hab2Diff = [];
Chrono302Diff = [];

meanFR30Same_Hab1 = [];
meanFR30Same_Trial1 = [];
meanFR30Same_Hab2 = [];
meanFR30Same_Trial2 = [];
meanFR30Same_Hab1Raw = [];
meanFR30Same_Trial1Raw = [];
meanFR30Same_Hab2Raw = [];
meanFR30Same_Trial2Raw = [];
meanFR30Same_Hab1_passed = [];
meanFR30Same_Trial1_passed = [];
meanFR30Same_Hab2_passed = [];
meanFR30Same_Trial2_passed = [];
meanFR30Same_Hab1All = [];
meanFR30Same_Trial1All = [];
meanFR30Same_Hab2All = [];
meanFR30Same_Trial2All = [];
meanFR30Same_Hab1All_passed = [];
meanFR30Same_Trial1All_passed = [];
meanFR30Same_Hab2All_passed = [];
meanFR30Same_Trial2All_passed = [];

meanFR30Diff_Hab1 = [];
meanFR30Diff_Trial1 = [];
meanFR30Diff_Hab2 = [];
meanFR30Diff_Trial2 = [];
meanFR30Diff_Hab1Raw = [];
meanFR30Diff_Trial1Raw = [];
meanFR30Diff_Hab2Raw = [];
meanFR30Diff_Trial2Raw = [];
meanFR30Diff_Hab1_passed = [];
meanFR30Diff_Trial1_passed = [];
meanFR30Diff_Hab2_passed = [];
meanFR30Diff_Trial2_passed = [];
meanFR30Diff_Hab1All = [];
meanFR30Diff_Trial1All = [];
meanFR30Diff_Hab2All = [];
meanFR30Diff_Trial2All = [];
meanFR30Diff_Hab1All_passed = [];
meanFR30Diff_Trial1All_passed = [];
meanFR30Diff_Hab2All_passed = [];
meanFR30Diff_Trial2All_passed = [];

meanFR120Same_Hab1 = [];
meanFR120Same_Trial1 = [];
meanFR120Same_Hab2 = [];
meanFR120Same_Trial2 = [];
meanFR120Same_Hab1Raw = [];
meanFR120Same_Trial1Raw = [];
meanFR120Same_Hab2Raw = [];
meanFR120Same_Trial2Raw = [];
meanFR120Same_Hab1_passed = [];
meanFR120Same_Trial1_passed = [];
meanFR120Same_Hab2_passed = [];
meanFR120Same_Trial2_passed = [];
meanFR120Same_Hab1All = [];
meanFR120Same_Trial1All = [];
meanFR120Same_Hab2All = [];
meanFR120Same_Trial2All = [];
meanFR120Same_Hab1All_passed = [];
meanFR120Same_Trial1All_passed = [];
meanFR120Same_Hab2All_passed = [];
meanFR120Same_Trial2All_passed = [];

meanFR120Diff_Hab1 = [];
meanFR120Diff_Trial1 = [];
meanFR120Diff_Hab2 = [];
meanFR120Diff_Trial2 = [];
meanFR120Diff_Hab1Raw = [];
meanFR120Diff_Trial1Raw = [];
meanFR120Diff_Hab2Raw = [];
meanFR120Diff_Trial2Raw = [];
meanFR120Diff_Hab1_passed = [];
meanFR120Diff_Trial1_passed = [];
meanFR120Diff_Hab2_passed = [];
meanFR120Diff_Trial2_passed = [];
meanFR120Diff_Hab1All = [];
meanFR120Diff_Trial1All = [];
meanFR120Diff_Hab2All = [];
meanFR120Diff_Trial2All = [];
meanFR120Diff_Hab1All_passed = [];
meanFR120Diff_Trial1All_passed = [];
meanFR120Diff_Hab2All_passed = [];
meanFR120Diff_Trial2All_passed = [];

bout1_2sec_30Same = [];
bout2_2sec_30Same = [];
boutHab1_2sec_30Same = [];
boutHab2_2sec_30Same = [];
bout1_2sec_30Diff = [];
bout2_2sec_30Diff = [];
boutHab1_2sec_30Diff = [];
boutHab2_2sec_30Diff = [];

bout1_3sec_30Same = [];
bout2_3sec_30Same = [];
boutHab1_3sec_30Same = [];
boutHab2_3sec_30Same = [];
bout1_3sec_30Diff = [];
bout2_3sec_30Diff = [];
boutHab1_3sec_30Diff = [];
boutHab2_3sec_30Diff = [];

bout1_32sec_30Same = [];
bout2_32sec_30Same = [];
boutHab1_32sec_30Same = [];
boutHab2_32sec_30Same = [];
bout1_32sec_30Diff = [];
bout2_32sec_30Diff = [];
boutHab1_32sec_30Diff = [];
boutHab2_32sec_30Diff = [];

bout1_2sec_120Same = [];
bout2_2sec_120Same = [];
boutHab1_2sec_120Same = [];
boutHab2_2sec_120Same = [];
bout1_2sec_120Diff = [];
bout2_2sec_120Diff = [];
boutHab1_2sec_120Diff = [];
boutHab2_2sec_120Diff = [];

bout1_3sec_120Same = [];
bout2_3sec_120Same = [];
boutHab1_3sec_120Same = [];
boutHab2_3sec_120Same = [];
bout1_3sec_120Diff = [];
bout2_3sec_120Diff = [];
boutHab1_3sec_120Diff = [];
boutHab2_3sec_120Diff = [];

bout1_32sec_120Same = [];
bout2_32sec_120Same = [];
boutHab1_32sec_120Same = [];
boutHab2_32sec_120Same = [];
bout1_32sec_120Diff = [];
bout2_32sec_120Diff = [];
boutHab1_32sec_120Diff = [];
boutHab2_32sec_120Diff = [];

bout1_2sec_30Same_CellReg = [];
bout2_2sec_30Same_CellReg = [];
boutHab1_2sec_30Same_CellReg = [];
boutHab2_2sec_30Same_CellReg = [];
bout1_2sec_30Diff_CellReg = [];
bout2_2sec_30Diff_CellReg = [];
boutHab1_2sec_30Diff_CellReg = [];
boutHab2_2sec_30Diff_CellReg = [];

bout1_2sec_120Same_CellReg = [];
bout2_2sec_120Same_CellReg = [];
boutHab1_2sec_120Same_CellReg = [];
boutHab2_2sec_120Same_CellReg = [];
bout1_2sec_120Diff_CellReg = [];
bout2_2sec_120Diff_CellReg = [];
boutHab1_2sec_120Diff_CellReg = [];
boutHab2_2sec_120Diff_CellReg = [];

boutHab1_2sec_30Diff_1h = [];
boutHab2_2sec_30Diff_1h = [];
boutHab1_2sec_30Same_1h = [];
boutHab2_2sec_30Same_1h = [];
boutHab1_2sec_30Diff_2h = [];
boutHab2_2sec_30Diff_2h = [];
boutHab1_2sec_30Same_2h = [];
boutHab2_2sec_30Same_2h = [];

boutHab1_2sec_120Diff_1h = [];
boutHab2_2sec_120Diff_1h = [];
boutHab1_2sec_120Same_1h = [];
boutHab2_2sec_120Same_1h = [];
boutHab1_2sec_120Diff_2h = [];
boutHab2_2sec_120Diff_2h = [];
boutHab1_2sec_120Same_2h = [];
boutHab2_2sec_120Same_2h = [];


for SITseshnum = 1 : length(SITproximityData.micenames)
    %     try
    mousename = SITproximityData.micenames{SITseshnum};
    folderparts = strsplit(SITproximityData.folderpaths{1,SITseshnum},'\');
    dateparts = strsplit(folderparts{end-3},' ');
    seshdate = dateparts{2};
    mkdir([cdold '\' seshdate '\' mousename])
    cd([SITproximityData.folderpaths{1,SITseshnum} '\..\..\ConcactenatedSession'])
    
    %Load in data
    load('ms.mat')
    load('frameMap1.mat')
    frameMap1 = frameMap;
    load('frameMap2.mat')
    frameMap2 = frameMap;
    load('SITstartFrame1.mat')
    startframe1 = startframe;
    ObjectStats1 = ObjectStats;
    load('SITstartFrame2.mat')
    startframe2 = startframe;
    ObjectStats2 = ObjectStats;
    load('badframes1.mat')
    bad1 = t;
    load('badframes2.mat')
    bad2 = t;
    load('HeadTrackingData1.mat')
    track1 = SINKdata;
    load('HeadTrackingData2.mat')
    track2 = SINKdata;
    
    cd([cdold '\' seshdate '\' mousename])
    
    %split sessions
    mins = 5;
    fps = 30;
    timeframes = mins*60*fps-1;
    
    track1 = track1(frameMap1,:);
    track2 = track2(frameMap2,:);
    trackhab1 = track1(bad1:timeframes+bad1,:);
    trackt1= track1(startframe1:timeframes+startframe1,:);
    trackhab2 = track2(bad2:timeframes+bad2,:);
    try
        trackt2= track2(startframe2:timeframes+startframe2,:);
    catch
        trackt2= track2(startframe2:end,:);
    end
    
    l1 = length(frameMap1);
    l2 = length(frameMap2);
    
    hab1 = ms.deconvolvedSig(bad1:timeframes+bad1,:)';
    trial1 = ms.deconvolvedSig(startframe1:timeframes+startframe1,:)';
    hab2 = ms.deconvolvedSig(l1+bad2:timeframes+bad2+l1,:)';
    try
        trial2 = ms.deconvolvedSig(startframe2+l1:timeframes+startframe2+l1,:)';
    catch
        trial2 = ms.deconvolvedSig(startframe2+l1:end,:)';
    end
    
    traces = cat(2,hab1,trial1,hab2,trial2);
    traces1 = cat(2,hab1,trial1);
    traces2 = cat(2,hab2,trial2);
    traces = traces(ms.exclude.SFPs,:);
    traces = normalize(traces,2,'range');
    downtraces = downsample(traces(:,9000:18000)',60);
    [maxval,maxframe]= max(downtraces,[],1);
    [~,newind] = sort(maxframe,'ascend');
    traces = traces(newind,:);
    
    figure
    imagesc(traces)
    vline(timeframes,'w:')
    vline(timeframes*2,'w:')
    vline(timeframes*3,'w:')
    
    for i = 1 : 2
        interactions = diff(SITproximityData.socialProxSessions{i,SITseshnum}.Object1Occ);
        interactionsHab = diff(SITproximityData.socialProxSessions{i,SITseshnum}.ObjectOccHab);
        
        enterzones = find(interactions > 0);
        exitzones = find(interactions < 0);
        enterzonesHab = find(interactionsHab > 0);
        exitzonesHab = find(interactionsHab < 0);
        
        if i == 1 && ~isempty(enterzonesHab) && ~isempty(exitzonesHab)
            mult = 0;
            vline(enterzonesHab,'g:')
            vline(exitzonesHab,'r:')
        elseif ~isempty(enterzonesHab) && ~isempty(exitzonesHab)
            mult = 1;
            vline(timeframes*(i)+enterzonesHab,'g:')
            vline(timeframes*(i)+exitzonesHab,'r:')
        end
        vline(timeframes*(i+mult)+enterzones,'g:')
        vline(timeframes*(i+mult)+exitzones,'r:')
        
        if i == 1
            enterzones1 = find(interactions > 0);
            exitzones1 = find(interactions < 0);
            enterzonesHab1 = find(interactionsHab > 0);
            exitzonesHab1 = find(interactionsHab < 0);
        else
            enterzones2 = find(interactions > 0);
            exitzones2 = find(interactions < 0);
            enterzonesHab2 = find(interactionsHab > 0);
            exitzonesHab2 = find(interactionsHab < 0);
        end
    end
    
    %% RateMap open field
    fieldsizex = 45;
    fieldsizey = 45;
    gridfactor = 2;
    celltotal = length(traces(:,1));
    map = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor),celltotal);
    map1 = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor),celltotal);
    map2 = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor),celltotal);
    mapHab1 = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor),celltotal);
    mapHab2 = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor),celltotal);
    mapframe = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    mapframe1 = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    mapframe2 = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    mapframeHab1 = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    mapframeHab2 = zeros(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    nanmask = ones(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    nanmask1 = ones(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    nanmask2 = ones(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    nanmaskHab1 = ones(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    nanmaskHab2 = ones(round(fieldsizey/gridfactor),round(fieldsizex/gridfactor));
    gcellind = find(ms.exclude.SFPs);
    
    track1 = cat(1,trackhab1,trackt1);
    track2 = cat(1,trackhab2,trackt2);
    
    track1 = track1./SITproximityData.socialProxSessions{1,SITseshnum}.pix2cm;
    track2 = track2./SITproximityData.socialProxSessions{2,SITseshnum}.pix2cm;
    
    trackhab1 = trackhab1./SITproximityData.socialProxSessions{1,SITseshnum}.pix2cm;
    trackhab2 = trackhab2./SITproximityData.socialProxSessions{2,SITseshnum}.pix2cm;
    
    trackt1 = trackt1./SITproximityData.socialProxSessions{1,SITseshnum}.pix2cm;
    trackt2 = trackt2./SITproximityData.socialProxSessions{2,SITseshnum}.pix2cm;
    
    min1x = min(track1(:,1));
    min1y = min(track1(:,2));
    
    min2x = min(track1(:,1));
    min2y = min(track1(:,2));
    
    track1(:,1) = track1(:,1)-min1x;
    track1(:,2) = track1(:,2)-min1y;
    trackhab1(:,1) = trackhab1(:,1)-min1x;
    trackhab1(:,2) = trackhab1(:,2)-min1y;
    trackt1(:,1) = trackt1(:,1)-min1x;
    trackt1(:,2) = trackt1(:,2)-min1y;
    
    track2(:,1) = track2(:,1)-min2x;
    track2(:,2) = track2(:,2)-min2y;
    trackhab2(:,1) = trackhab2(:,1)-min2x;
    trackhab2(:,2) = trackhab2(:,2)-min2y;
    trackt2(:,1) = trackt2(:,1)-min2x;
    trackt2(:,2) = trackt2(:,2)-min2y;
    
    for i = gridfactor : gridfactor : fieldsizex
        for j = gridfactor : gridfactor : fieldsizey
            if i == gridfactor && (j ~= fieldsizey && j ~= gridfactor)
                t1 = find(trackt1(:,1) < i & trackt1(:,2) < j & trackt1(:,2) >= j-gridfactor);
                t2 = find(trackt2(:,1) < i & trackt2(:,2) < j & trackt2(:,2) >= j-gridfactor);
                habt1 = find(trackhab1(:,1) < i & trackhab1(:,2) < j & trackhab1(:,2) >= j-gridfactor);
                habt2 = find(trackhab2(:,1) < i & trackhab2(:,2) < j & trackhab2(:,2) >= j-gridfactor);
            elseif i == fieldsizex && (j ~= fieldsizey && j ~= gridfactor)
                t1 = find(trackt1(:,1) >= i-gridfactor & trackt1(:,2) < j & trackt1(:,2) >= j-gridfactor);
                t2 = find(trackt2(:,1) >= i-gridfactor & trackt2(:,2) < j & trackt2(:,2) >= j-gridfactor);
                habt1 = find(trackhab1(:,1) >= i-gridfactor & trackhab1(:,2) < j & trackhab1(:,2) >= j-gridfactor);
                habt2 = find(trackhab2(:,1) >= i-gridfactor & trackhab2(:,2) < j & trackhab2(:,2) >= j-gridfactor);
            elseif j == gridfactor && (i ~= fieldsizex && i ~= gridfactor)
                t1 = find(trackt1(:,2) < j & trackt1(:,1) < i & trackt1(:,1) >= i-gridfactor);
                t2 = find(trackt2(:,2) < j & trackt2(:,1) < i & trackt2(:,1) >= i-gridfactor);
                habt1 = find(trackhab1(:,2) < j & trackhab1(:,1) < i & trackhab1(:,1) >= i-gridfactor);
                habt2 = find(trackhab2(:,2) < j & trackhab2(:,1) < i & trackhab2(:,1) >= i-gridfactor);
            elseif j == fieldsizey && (i ~= fieldsizex && i ~= gridfactor)
                t1 = find(trackt1(:,2) >= j-gridfactor & trackt1(:,1) < i & trackt1(:,1) >= i-gridfactor);
                t2 = find(trackt2(:,2) >= j-gridfactor & trackt2(:,1) < i & trackt2(:,1) >= i-gridfactor);
                habt1 = find(trackhab1(:,2) >= j-gridfactor & trackhab1(:,1) < i & trackhab1(:,1) >= i-gridfactor);
                habt2 = find(trackhab2(:,2) >= j-gridfactor & trackhab2(:,1) < i & trackhab2(:,1) >= i-gridfactor);
            elseif i == gridfactor && j == gridfactor
                t1 = find(trackt1(:,2) < j & trackt1(:,1) < i);
                t2 = find(trackt2(:,2) < j & trackt2(:,1) < i);
                habt1 = find(trackhab1(:,2) < j & trackhab1(:,1) < i);
                habt2 = find(trackhab2(:,2) < j & trackhab2(:,1) < i);
            elseif i == fieldsizex && j == fieldsizey
                t1 = find(trackt1(:,2) >= j-gridfactor & trackt1(:,1) >= i-gridfactor);
                t2 = find(trackt2(:,2) >= j-gridfactor & trackt2(:,1) >= i-gridfactor);
                habt1 = find(trackhab1(:,2) >= j-gridfactor & trackhab1(:,1) >= i-gridfactor);
                habt2 = find(trackhab2(:,2) >= j-gridfactor & trackhab2(:,1) >= i-gridfactor);
            elseif i == gridfactor && j == fieldsizey
                t1 = find(trackt1(:,1) < i & trackt1(:,2) >= j-gridfactor);
                t2 = find(trackt2(:,1) < i & trackt2(:,2) >= j-gridfactor);
                habt1 = find(trackhab1(:,1) < i & trackhab1(:,2) >= j-gridfactor);
                habt2 = find(trackhab2(:,1) < i & trackhab2(:,2) >= j-gridfactor);
            elseif j == gridfactor && i == fieldsizex
                t1 = find(trackt1(:,2) < j & trackt1(:,1) >= i-gridfactor);
                t2 = find(trackt2(:,2) < j & trackt2(:,1) >= i-gridfactor);
                habt1 = find(trackhab1(:,2) < j & trackhab1(:,1) >= i-gridfactor);
                habt2 = find(trackhab2(:,2) < j & trackhab2(:,1) >= i-gridfactor);
            else
                t1 = intersect(find(trackt1(:,1) < i & trackt1(:,1) >= i-gridfactor), find(trackt1(:,2) < j & trackt1(:,2) >= j-gridfactor));
                t2 = intersect(find(trackt2(:,1) < i & trackt2(:,1) >= i-gridfactor), find(trackt2(:,2) < j & trackt2(:,2) >= j-gridfactor));
                habt1 = intersect(find(trackhab1(:,1) < i & trackhab1(:,1) >= i-gridfactor), find(trackhab1(:,2) < j & trackhab1(:,2) >= j-gridfactor));
                habt2 = intersect(find(trackhab2(:,1) < i & trackhab2(:,1) >= i-gridfactor), find(trackhab2(:,2) < j & trackhab2(:,2) >= j-gridfactor));
            end
            if ~(isempty(t1) && isempty(t2) && isempty(habt1) && isempty(habt2))
                for cellnum = 1 : celltotal
                    map(round(i/gridfactor),round(j/gridfactor),cellnum) = sum(trial1(gcellind(cellnum),t1));
                    map(round(i/gridfactor),round(j/gridfactor),cellnum) = sum(hab1(gcellind(cellnum),habt1)) + sum(trial2(gcellind(cellnum),t2)) + sum(hab2(gcellind(cellnum),habt2)) + map(round(i/gridfactor),round(j/gridfactor),cellnum);
                    if ~isempty(t1)
                        map1(round(i/gridfactor),round(j/gridfactor),cellnum) = sum(trial1(gcellind(cellnum),t1));
                    end
                    if ~isempty(t2)
                        map2(round(i/gridfactor),round(j/gridfactor),cellnum) = sum(trial2(gcellind(cellnum),t2));
                    end
                    if ~isempty(habt1)
                        mapHab1(round(i/gridfactor),round(j/gridfactor),cellnum) = sum(hab1(gcellind(cellnum),habt1));
                    end
                    if ~isempty(habt2)
                        mapHab2(round(i/gridfactor),round(j/gridfactor),cellnum) = sum(hab2(gcellind(cellnum),habt2));
                    end
                end
                if ~isempty(t1)
                    mapframe(round(i/gridfactor),round(j/gridfactor)) = length(t1);
                    mapframe1(round(i/gridfactor),round(j/gridfactor)) = length(t1);
                else
                    nanmask1(round(i/gridfactor),round(j/gridfactor)) = nan;
                end
                if ~isempty(t2)
                    mapframe(round(i/gridfactor),round(j/gridfactor)) = length(t2);
                    mapframe2(round(i/gridfactor),round(j/gridfactor)) = length(t2);
                else
                    nanmask2(round(i/gridfactor),round(j/gridfactor)) = nan;
                end
                if ~isempty(habt1)
                    mapframe(round(i/gridfactor),round(j/gridfactor)) = length(habt1);
                    mapframeHab1(round(i/gridfactor),round(j/gridfactor)) = length(habt1);
                else
                    nanmaskHab1(round(i/gridfactor),round(j/gridfactor)) = nan;
                end
                if ~isempty(habt2)
                    mapframe(round(i/gridfactor),round(j/gridfactor)) = length(habt2);
                    mapframeHab2(round(i/gridfactor),round(j/gridfactor)) = length(habt2);
                else
                    nanmaskHab2(round(i/gridfactor),round(j/gridfactor)) = nan;
                end
            else
                nanmask(round(i/gridfactor),round(j/gridfactor)) = nan;
                nanmask1(round(i/gridfactor),round(j/gridfactor)) = nan;
                nanmask2(round(i/gridfactor),round(j/gridfactor)) = nan;
                nanmaskHab1(round(i/gridfactor),round(j/gridfactor)) = nan;
                nanmaskHab2(round(i/gridfactor),round(j/gridfactor)) = nan;
            end
        end
    end
    
    map = map./mapframe;
    map1 = map1./mapframe1;
    map2 = map2./mapframe2;
    mapHab1 = mapHab1./mapframeHab1;
    mapHab2 = mapHab2./mapframeHab2;
    
    map = rot90(map);
    nanmask = rot90(nanmask);
    map1 = rot90(map1);
    nanmask1 = rot90(nanmask1);
    map2 = rot90(map2);
    nanmask2 = rot90(nanmask2);
    mapframe1 = rot90(mapframe1);
    mapframe2 = rot90(mapframe2);
    
    mapHab1 = rot90(mapHab1);
    nanmaskHab1 = rot90(nanmaskHab1);
    mapHab2 = rot90(mapHab2);
    nanmaskHab2 = rot90(nanmaskHab2);
    mapframeHab1 = rot90(mapframeHab1);
    mapframeHab2 = rot90(mapframeHab2);
    
    filtkernel = 1;
    count = 1;
    %{
        figure('Position', get(0, 'Screensize'));
        for i = 1 : celltotal
            map1(isnan(map1)) = 0;
            map2(isnan(map2)) = 0;
            map1(:,:,i) = imgaussfilt(map1(:,:,i),filtkernel).*nanmask1;
            map2(:,:,i) = imgaussfilt(map2(:,:,i),filtkernel).*nanmask2;
            
            mapHab1(isnan(mapHab1)) = 0;
            mapHab2(isnan(mapHab2)) = 0;
            mapHab1(:,:,i) = imgaussfilt(mapHab1(:,:,i),filtkernel).*nanmaskHab1;
            mapHab2(:,:,i) = imgaussfilt(mapHab2(:,:,i),filtkernel).*nanmaskHab2;
            
            if mod(count,2) == 0 && count ~=1
                count = count + 1;
            end
            subplot(4,4,count)
            imagesc(mapHab1(:,:,i),'AlphaData',~isnan(nanmaskHab1))
            colormap('turbo')
            caxis([0 max(nanmax(mapHab1(:,:,i)))])
            title(['Habituation 1 Cell ' num2str(gcellind(i))])
            subplot(4,4,count+4)
            plot(trackhab1(:,1),trackhab1(:,2))
            hold on
            indfiring = find(hab1(gcellind(i),:));
            scatter(trackhab1(indfiring,1),trackhab1(indfiring,2),'r.')
            set(gca,'visible','off')
            
            subplot(4,4,count+1)
            imagesc(map1(:,:,i),'AlphaData',~isnan(nanmask1))
            colormap('turbo')
            caxis([0 max(nanmax(map1(:,:,i)))])
            title(['Trial1 Cell ' num2str(gcellind(i))])
            subplot(4,4,count+5)
            plot(trackt1(:,1),trackt1(:,2))
            hold on
            indfiring = find(trial1(gcellind(i),:));
            scatter(trackt1(indfiring,1),trackt1(indfiring,2),'r.')
            set(gca,'visible','off')
            
            subplot(4,4,count+2)
            imagesc(mapHab2(:,:,i),'AlphaData',~isnan(nanmaskHab2))
            colormap('turbo')
            caxis([0 max(nanmax(mapHab2(:,:,i)))])
            title(['Habituation 2 Cell ' num2str(gcellind(i))])
            subplot(4,4,count+6)
            plot(trackhab2(:,1),trackhab2(:,2))
            hold on
            indfiring = find(hab2(gcellind(i),:));
            scatter(trackhab2(indfiring,1),trackhab2(indfiring,2),'r.')
            set(gca,'visible','off')
            
            subplot(4,4,count+3)
            imagesc(map2(:,:,i),'AlphaData',~isnan(nanmask2))
            colormap('turbo')
            caxis([0 max(nanmax(map2(:,:,i)))])
            title(['Trial2 Cell ' num2str(gcellind(i))])
            subplot(4,4,count+7)
            plot(trackt2(:,1),trackt2(:,2))
            hold on
            indfiring = find(trial2(gcellind(i),:));
            scatter(trackt2(indfiring,1),trackt2(indfiring,2),'r.')
            set(gca,'visible','off')
            sgtitle([seshdate ' ' mousename ' ' SITnovelty_mins{SITseshnum,1} ' ' num2str(SITnovelty_mins{SITseshnum,2}) 'mins'])
            if count == 1
                count = 9;
            else
                count = count+1;
            end
            if count > 9
                count = 1;
                saveas(gcf,['FiringRateMap_Trial1vsTrial2_' num2str(gcellind(i))])
                saveas(gcf,['FiringRateMap_Trial1vsTrial2_' num2str(gcellind(i)) '.jpeg'])
                pause(0.01)
                clf
            end
        end
    %}
    
    %% Mean Interaction Activity
    % ztrial = traces
    ztrials = normalize(cat(2,hab1,trial1,hab2,trial2),2);
    zhab1 = ztrials(:,1:length(hab1(1,:)));
    ztrial1 = ztrials(:,length(zhab1(1,:))+1:length(zhab1(1,:))+length(trial1(1,:)));
    zhab2 = ztrials(:,length(traces1(1,:))+1:length(traces1(1,:))+length(hab2(1,:)));
    ztrial2 = ztrials(:,length(traces1(1,:))+length(hab2(1,:))+1:end);
    ztrials_2 = normalize(cat(2,hab1(:,4501:end),trial1,hab2(:,4501:end),trial2),2);
    
    if ~isempty(exitzones1) && exitzones1(1) < enterzones1(1)
        exitzones1(1) = [];
    end
    if ~isempty(exitzones2) && exitzones2(1) < enterzones2(1)
        exitzones2(1) = [];
    end
    if ~isempty(exitzonesHab1) && exitzonesHab1(1) < enterzonesHab1(1)
        exitzonesHab1(1) = [];
    end
    if ~isempty(exitzonesHab2) && exitzonesHab2(1) < enterzonesHab2(1)
        exitzonesHab2(1) = [];
    end
    
    if length(enterzones1) ~= length(exitzones1)
        if length(enterzones1) > length(exitzones1)
            enterzones1(end) = [];
        else
            if exitzones1(1) < enterzones1(1)
                exitzones1(1) = [];
            else
                exitzones1(end) = [];
            end
        end
    end
    if length(enterzones2) ~= length(exitzones2)
        if length(enterzones2) > length(exitzones2)
            enterzones2(end) = [];
        else
            if exitzones2(1) < enterzones2(1)
                exitzones2(1) = [];
            else
                exitzones2(end) = [];
            end
        end
    end
    if length(enterzonesHab1) ~= length(exitzonesHab1)
        if length(enterzonesHab1) > length(exitzonesHab1)
            enterzonesHab1(end) = [];
        else
            if exitzonesHab1(1) < enterzonesHab1(1)
                exitzonesHab1(1) = [];
            else
                exitzonesHab1(end) = [];
            end
        end
    end
    if length(enterzonesHab2) ~= length(exitzonesHab2)
        if length(enterzonesHab2) > length(exitzonesHab2)
            enterzonesHab2(end) = [];
        else
            if exitzonesHab2(1) < enterzonesHab2(1)
                exitzonesHab2(1) = [];
            else
                exitzonesHab2(end) = [];
            end
        end
    end        
    
    interactionlength1 = cat(1,interactionlength1, exitzones1 - enterzones1);
    interactionlength2 = cat(1,interactionlength2, exitzones2 - enterzones2);
    interactionlengthHab1 = cat(1,interactionlengthHab1, exitzonesHab1 - enterzonesHab1);
    interactionlengthHab2 = cat(1,interactionlengthHab2, exitzonesHab2 - enterzonesHab2);
    
    enterzonesall{SITseshnum,1} = enterzones1;
    enterzonesall{SITseshnum,2} = enterzones2;
    enterzonesall{SITseshnum,3} = enterzonesHab1;
    enterzonesall{SITseshnum,4} = enterzonesHab2;       
    
    exitzonesall{SITseshnum,1} = exitzones1;
    exitzonesall{SITseshnum,2} = exitzones2;
    exitzonesall{SITseshnum,3} = exitzonesHab1;
    exitzonesall{SITseshnum,4} = exitzonesHab2;
    
    timebinnum = 5;
    interactionsMean = zeros(celltotal,timebinnum);
    zoneActivitysplit1 = zeros(celltotal,length(enterzones1),timebinnum);
    outzoneActivitysplit1 = zeros(celltotal);
    zoneActivityRaw1 = zeros(celltotal,length(enterzones1));
    zoneActivityRaw2 = zeros(celltotal,length(enterzones2));
    zoneActivityHabRaw1 = zeros(celltotal,length(enterzonesHab1));
    zoneActivityHabRaw2 = zeros(celltotal,length(enterzonesHab2));
    
    zoneActivity1 = zeros(celltotal,length(enterzones1));
    zoneActivity2 = zeros(celltotal,length(enterzones2));
    zoneActivityHab1 = zeros(celltotal,length(enterzonesHab1));
    zoneActivityHab2 = zeros(celltotal,length(enterzonesHab2));        
    
    indshortint = [];
    
    for time = 1 : timebinnum
        for i = 1 : length(enterzones1)
            if (exitzones1(i)-enterzones1(i)) >= 60
                binsize = round((exitzones1(i)-enterzones1(i))/timebinnum);
                if binsize > 1
                    for j = 1 : celltotal
                        zoneActivitysplit1(j,i,time) = sum(ztrial1(gcellind(j),enterzones1(i)+binsize*(time-1):enterzones1(i)+binsize*(time)))/binsize;
                        zoneActivityRaw1(j,i) = sum(trial1(gcellind(j),enterzones1(i):exitzones1(i)))/(exitzones1(i)-enterzones1(i));
                        zoneActivity1(j,i) = sum(ztrial1(gcellind(j),enterzones1(i):exitzones1(i)))/(exitzones1(i)-enterzones1(i));
                    end
                end
            else
                indshortint = cat(1,indshortint,i);
            end
        end
    end
    indshortint = sort(unique(indshortint),'descend');
    zoneActivityRaw1(:,indshortint) = [];
    zoneActivity1(:,indshortint) = [];
    
    indshortint = [];
    indkeep = [];
    hab1_2ind = [];
    zoneActivitysplitHab1 = zeros(celltotal,length(enterzonesHab1),timebinnum);
    for time = 1 : timebinnum
        for i = 1 : length(enterzonesHab1)
            if (exitzonesHab1(i)-enterzonesHab1(i)) >= 60
                binsize = round((exitzonesHab1(i)-enterzonesHab1(i))/timebinnum);
                if binsize > 1
                    for j = 1 : celltotal
                        zoneActivitysplitHab1(j,i,time) = sum(zhab1(gcellind(j),enterzonesHab1(i)+binsize*(time-1):enterzonesHab1(i)+binsize*(time)))/binsize;
                        zoneActivityHabRaw1(j,i) = sum(hab1(gcellind(j),enterzonesHab1(i):exitzonesHab1(i)))/(exitzonesHab1(i)-enterzonesHab1(i));
                        zoneActivityHab1(j,i) = sum(zhab1(gcellind(j),enterzonesHab1(i):exitzonesHab1(i)))/(exitzonesHab1(i)-enterzonesHab1(i));                        
                    end
                    indkeep = cat(1,indkeep,i);
                end
            else
                indshortint = cat(1,indshortint,i);
            end
            if enterzonesHab1(i) <= 4500
                hab1_2ind = i;
            end
        end
    end
    indshortint = sort(unique(indshortint),'descend');
    zoneActivityHabRaw1(:,indshortint) = [];
    zoneActivityHab1(:,indshortint) = [];
    
    indkeep = unique(indkeep);
    
    if isempty(hab1_2ind)        
        hab1_2ind = 0;
        enterzonesHab1_2 = enterzonesHab1-4501;
        exitzonesHab1_2 = exitzonesHab1-4501;
    else        
        enterzonesHab1_2 = enterzonesHab1(hab1_2ind+1:end)-4501;
        exitzonesHab1_2 = exitzonesHab1(hab1_2ind+1:end)-4501;        
        if isempty(find(indkeep<=hab1_2ind))
            hab1_2ind = 0;
        else
            hab1_2ind = max(find(indkeep<=hab1_2ind));
        end
    end                    
        
    zoneActivityHab1_2 = zoneActivityHab1(:,hab1_2ind+1:end);
    
    indshortint = [];    
    zoneActivitysplit2 = zeros(celltotal,length(enterzones2),timebinnum);    
    for time = 1 : timebinnum
        for i = 1 : length(enterzones2)
            if (exitzones2(i)-enterzones2(i)) >= 60
            binsize = round((exitzones2(i)-enterzones2(i))/timebinnum);
            if binsize > 1
                for j = 1 : celltotal
                    zoneActivitysplit2(j,i,time) = sum(ztrial2(gcellind(j),enterzones2(i)+binsize*(time-1):enterzones2(i)+binsize*(time)))/binsize;
                    zoneActivityRaw2(j,i) = sum(trial2(gcellind(j),enterzones2(i):exitzones2(i)))/(exitzones2(i)-enterzones2(i));
                    zoneActivity2(j,i) = sum(ztrial2(gcellind(j),enterzones2(i):exitzones2(i)))/(exitzones2(i)-enterzones2(i));                    
                end
            end            
            else
                indshortint = cat(1,indshortint,i);
            end            
        end
    end
    indshortint = sort(unique(indshortint),'descend');
    zoneActivityRaw2(:,indshortint) = [];
    zoneActivity2(:,indshortint) = [];      
    
    indshortint = [];
    
    indkeep = [];
    hab2_2ind = [];
    zoneActivitysplitHab2 = zeros(celltotal,length(enterzonesHab2),timebinnum);
    for time = 1 : timebinnum
        for i = 1 : length(enterzonesHab2)
            if (exitzonesHab2(i)-enterzonesHab2(i)) >= 60
            binsize = round((exitzonesHab2(i)-enterzonesHab2(i))/timebinnum);
            if binsize > 1
                for j = 1 : celltotal
                    zoneActivitysplitHab2(j,i,time) = sum(zhab2(gcellind(j),enterzonesHab2(i)+binsize*(time-1):enterzonesHab2(i)+binsize*(time)))/binsize;
                    zoneActivityHabRaw2(j,i) = sum(hab2(gcellind(j),enterzonesHab2(i):exitzonesHab2(i)))/(exitzonesHab2(i)-enterzonesHab2(i));
                    zoneActivityHab2(j,i) = sum(zhab2(gcellind(j),enterzonesHab2(i):exitzonesHab2(i)))/(exitzonesHab2(i)-enterzonesHab2(i));
                end                
                indkeep = cat(1,indkeep,i);
            end
            else
                indshortint = cat(1,indshortint,i);
            end
            if enterzonesHab2(i) <= 4500
                    hab2_2ind = i;
            end
        end
    end
    indshortint = sort(unique(indshortint),'descend');
    zoneActivityHabRaw2(:,indshortint) = [];
    zoneActivityHab2(:,indshortint) = [];
    
    indkeep = unique(indkeep);
    if isempty(hab2_2ind)        
        hab2_2ind = 0;
        enterzonesHab2_2 = enterzonesHab2-4501;
        exitzonesHab2_2 = exitzonesHab2-4501;
    else
        enterzonesHab2_2 = enterzonesHab2(hab2_2ind+1:end)-4501;
        exitzonesHab2_2 = exitzonesHab2(hab2_2ind+1:end)-4501;
        if isempty(find(indkeep<=hab2_2ind))
            hab2_2ind = 0;
        else
            hab2_2ind = max(find(indkeep<=hab2_2ind));
        end
    end         
        
    zoneActivityHab2_2 = zoneActivityHab2(:,hab2_2ind+1:end);
    
%     if isempty(zoneActivityHab1_2)
% %         zoneActivityHab1_2 = zeros(length(zoneActivityHab1(:,1)),1);
%     end
%     if isempty(zoneActivityHab2_2)
% %         zoneActivityHab2_2 = zeros(length(zoneActivityHab2(:,1)),1);
%     end
    
    timebin2 = 50;
    ChronoActivity1 = zeros(celltotal,timebin2);
    ChronoActivity2 = zeros(celltotal,timebin2);
    ChronoActivityHab1 = zeros(celltotal,timebin2);
    ChronoActivityHab2 = zeros(celltotal,timebin2);        
    
    catActivity1 = [];
    catActivity1Raw = [];
    for interact = 1 : length(enterzones1)
        if exitzones1(interact) - enterzones1(interact) >= secframes2            
            if interact == 1
                catActivity1 = ztrial1(gcellind,enterzones1(interact):exitzones1(interact));
                catActivity1Raw = trial1(gcellind,enterzones1(interact):exitzones1(interact));
            else
                catActivity1 = cat(2,catActivity1,ztrial1(gcellind,enterzones1(interact):exitzones1(interact)));
                catActivity1Raw = cat(2,catActivity1Raw,trial1(gcellind,enterzones1(interact):exitzones1(interact)));
            end
        end
    end
    
    catActivity2 = [];
    catActivity2Raw = [];
    for interact = 1 : length(enterzones2)
        if exitzones2(interact) - enterzones2(interact) >= secframes2
            if interact == 1
                catActivity2 = ztrial2(gcellind,enterzones1(interact):exitzones1(interact));
                catActivity2Raw = trial2(gcellind,enterzones1(interact):exitzones1(interact));
            else
                catActivity2 = cat(2,catActivity2,ztrial2(gcellind,enterzones2(interact):exitzones2(interact)));
                catActivity2Raw = cat(2,catActivity2Raw,trial2(gcellind,enterzones2(interact):exitzones2(interact)));
            end
        end
    end
    
    catActivityHab1 = [];
    catActivityHab1Raw = [];
    for interact = 1 : length(enterzonesHab1)
        if exitzonesHab1(interact) - enterzonesHab1(interact) >= secframes2            
            if interact == 1
                catActivityHab1 = zhab1(gcellind,enterzones1(interact):exitzonesHab1(interact));
                catActivityHab1Raw = hab1(gcellind,enterzones1(interact):exitzonesHab1(interact));
            else
                catActivityHab1 = cat(2,catActivityHab1,zhab1(gcellind,enterzonesHab1(interact):exitzonesHab1(interact)));
                catActivityHab1Raw = cat(2,catActivityHab1Raw,hab1(gcellind,enterzonesHab1(interact):exitzonesHab1(interact)));
            end
        end
    end
    
    catActivityHab2 = [];
    catActivityHab2Raw = [];
    for interact = 1 : length(enterzonesHab2)
        if exitzonesHab2(interact) - enterzonesHab2(interact) >= secframes2            
            if interact == 1
                catActivityHab2 = zhab2(gcellind,enterzonesHab2(interact):exitzonesHab2(interact));
                catActivityHab2Raw = hab2(gcellind,enterzonesHab2(interact):exitzonesHab2(interact));
            else
                catActivityHab2 = cat(2,catActivityHab2,zhab2(gcellind,enterzonesHab2(interact):exitzonesHab2(interact)));
                catActivityHab2Raw = cat(2,catActivityHab2Raw,hab2(gcellind,enterzonesHab2(interact):exitzonesHab2(interact)));
            end
        end
    end
    if length(catActivity1) == 0
        catActivity1 = zeros(celltotal,timebin2);
        catActivity1Raw = catActivity1;
    end
    if length(catActivity2) == 0
        catActivity2 = zeros(celltotal,timebin2);
        catActivity2Raw = catActivity2;
    end
    if length(catActivityHab1) == 0
        catActivityHab1 = zeros(celltotal,timebin2);
        catActivityHab1Raw = catActivityHab1;
    end
    if length(catActivityHab2) == 0
        catActivityHab2 = zeros(celltotal,timebin2);
        catActivityHab2Raw = catActivityHab2;
    end
    
    %bout analysis
    secframes2 = 60;
    secframes3 = 90;
    
    %bout1
    for i = 1: length(enterzones1)
        if enterzones1(i) > secframes2
            if exitzones1(i) - enterzones1(i) >= secframes2
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes2
                            bout1_2sec_30Same = cat(1,bout1_2sec_30Same,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes2));
                        elseif i == 1
                            bout1_2sec_30Same = cat(1,bout1_2sec_30Same,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes2));
                        end
%                         bout1_2sec_30Same_1h = boutcat_1h(enterzones1,exitzones1, secframes2,gcellind,ztrial1,bout1_2sec_30Same_1h,i);
%                         bout1_2sec_30Same_2h = boutcat_2h(enterzones1,exitzones1, secframes2,gcellind,ztrial1,bout1_2sec_30Same_2h,i);
                    else
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes2
                            bout1_2sec_120Same = cat(1,bout1_2sec_120Same,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes2));
                        elseif i == 1
                            bout1_2sec_120Same = cat(1,bout1_2sec_120Same,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes2));
                        end
%                         bout1_2sec_120Same_1h = boutcat_1h(enterzones1,exitzones1, secframes2,gcellind,ztrial1,bout1_2sec_120Same_1h,i);
%                         bout1_2sec_120Same_2h = boutcat_2h(enterzones1,exitzones1, secframes2,gcellind,ztrial1,bout1_2sec_120Same_2h,i);
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes2
                            bout1_2sec_30Diff = cat(1,bout1_2sec_30Diff,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes2));
                        elseif i == 1
                            bout1_2sec_30Diff = cat(1,bout1_2sec_30Diff,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes2));
                        end
%                         bout1_2sec_30Diff_1h = boutcat_1h(enterzones1,exitzones1, secframes2,gcellind,ztrial1,bout1_2sec_30Diff_1h,i);
%                         bout1_2sec_30Diff_2h = boutcat_2h(enterzones1,exitzones1, secframes2,gcellind,ztrial1,bout1_2sec_30Diff_2h,i);
                    else
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes2
                            bout1_2sec_120Diff = cat(1,bout1_2sec_120Diff,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes2));
                        elseif i == 1
                            bout1_2sec_120Diff = cat(1,bout1_2sec_120Diff,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes2));
                        end
%                         bout1_2sec_120Diff_1h = boutcat_1h(enterzones1,exitzones1, secframes2,gcellind,ztrial1,bout1_2sec_120Diff_1h,i);
%                         bout1_2sec_120Diff_2h = boutcat_2h(enterzones1,exitzones1, secframes2,gcellind,ztrial1,bout1_2sec_120Diff_2h,i);
                    end
                end
            end
        end
        if enterzones1(i) > secframes3
            if exitzones1(i) - enterzones1(i) >= secframes3
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes3
                            bout1_3sec_30Same = cat(1,bout1_3sec_30Same,ztrial1(gcellind,enterzones1(i)-secframes3:enterzones1(i)+secframes3));
                        elseif i == 1
                            bout1_3sec_30Same = cat(1,bout1_3sec_30Same,ztrial1(gcellind,enterzones1(i)-secframes3:enterzones1(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes3
                            bout1_3sec_120Same = cat(1,bout1_3sec_120Same,ztrial1(gcellind,enterzones1(i)-secframes3:enterzones1(i)+secframes3));
                        elseif i == 1
                            bout1_3sec_120Same = cat(1,bout1_3sec_120Same,ztrial1(gcellind,enterzones1(i)-secframes3:enterzones1(i)+secframes3));
                        end
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes3
                            bout1_3sec_30Diff = cat(1,bout1_3sec_30Diff,ztrial1(gcellind,enterzones1(i)-secframes3:enterzones1(i)+secframes3));
                        elseif i == 1
                            bout1_3sec_30Diff = cat(1,bout1_3sec_30Diff,ztrial1(gcellind,enterzones1(i)-secframes3:enterzones1(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes3
                            bout1_3sec_120Diff = cat(1,bout1_3sec_120Diff,ztrial1(gcellind,enterzones1(i)-secframes3:enterzones1(i)+secframes3));
                        elseif i == 1
                            bout1_3sec_120Diff = cat(1,bout1_3sec_120Diff,ztrial1(gcellind,enterzones1(i)-secframes3:enterzones1(i)+secframes3));
                        end
                    end
                end
                % 2 pre and 3 post interactions
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes2
                            bout1_32sec_30Same = cat(1,bout1_32sec_30Same,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes3));
                        elseif i == 1
                            bout1_32sec_30Same = cat(1,bout1_32sec_30Same,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes2
                            bout1_32sec_120Same = cat(1,bout1_32sec_120Same,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes3));
                        elseif i == 1
                            bout1_32sec_120Same = cat(1,bout1_32sec_120Same,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes3));
                        end
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes2
                            bout1_32sec_30Diff = cat(1,bout1_32sec_30Diff,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes3));
                        elseif i == 1
                            bout1_32sec_30Diff = cat(1,bout1_32sec_30Diff,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzones1(i) - exitzones1(i-1) >= secframes2
                            bout1_32sec_120Diff = cat(1,bout1_32sec_120Diff,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes3));
                        elseif i == 1
                            bout1_32sec_120Diff = cat(1,bout1_32sec_120Diff,ztrial1(gcellind,enterzones1(i)-secframes2:enterzones1(i)+secframes3));
                        end
                    end
                end
            end
        end
    end
    
    % bout2
    for i = 1: length(enterzones2)
        if enterzones2(i) > secframes2
            if exitzones2(i) - enterzones2(i) >= secframes2
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes2
                            bout2_2sec_30Same = cat(1,bout2_2sec_30Same,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes2));
                        elseif i == 1
                            bout2_2sec_30Same = cat(1,bout2_2sec_30Same,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes2));
                        end
%                         bout2_2sec_30Same_1h = boutcat_1h(enterzones2,exitzones2, secframes2,gcellind,ztrial2,bout2_2sec_30Same_1h,i);
%                         bout2_2sec_30Same_2h = boutcat_2h(enterzones2,exitzones2, secframes2,gcellind,ztrial2,bout2_2sec_30Same_2h,i);
                    else
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes2
                            bout2_2sec_120Same = cat(1,bout2_2sec_120Same,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes2));
                        elseif i == 1
                            bout2_2sec_120Same = cat(1,bout2_2sec_120Same,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes2));
                        end
%                         bout2_2sec_120Same_1h = boutcat_1h(enterzones2,exitzones2, secframes2,gcellind,ztrial2,bout2_2sec_120Same_1h,i);
%                         bout2_2sec_120Same_2h = boutcat_2h(enterzones2,exitzones2, secframes2,gcellind,ztrial2,bout2_2sec_120Same_2h,i);
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes2
                            bout2_2sec_30Diff = cat(1,bout2_2sec_30Diff,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes2));
                        elseif i == 1
                            bout2_2sec_30Diff = cat(1,bout2_2sec_30Diff,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes2));
                        end
%                         bout2_2sec_30Diff_1h = boutcat_1h(enterzones2,exitzones2, secframes2,gcellind,ztrial2,bout2_2sec_30Diff_1h,i);
%                         bout2_2sec_30Diff_2h = boutcat_2h(enterzones2,exitzones2, secframes2,gcellind,ztrial2,bout2_2sec_30Diff_2h,i);
                    else
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes2
                            bout2_2sec_120Diff = cat(1,bout2_2sec_120Diff,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes2));
                        elseif i == 1
                            bout2_2sec_120Diff = cat(1,bout2_2sec_120Diff,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes2));
                        end
%                         bout2_2sec_120Diff_1h = boutcat_1h(enterzones2,exitzones2, secframes2,gcellind,ztrial2,bout2_2sec_120Diff_1h,i);
%                         bout2_2sec_120Diff_2h = boutcat_2h(enterzones2,exitzones2, secframes2,gcellind,ztrial2,bout2_2sec_120Diff_2h,i);
                    end
                end
            end
        end
        if enterzones2(i) > secframes3
            if exitzones2(i) - enterzones2(i) >= secframes3
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes3
                            bout2_3sec_30Same = cat(1,bout2_3sec_30Same,ztrial2(gcellind,enterzones2(i)-secframes3:enterzones2(i)+secframes3));
                        elseif i == 1
                            bout2_3sec_30Same = cat(1,bout2_3sec_30Same,ztrial2(gcellind,enterzones2(i)-secframes3:enterzones2(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes3
                            bout2_3sec_120Same = cat(1,bout2_3sec_120Same,ztrial2(gcellind,enterzones2(i)-secframes3:enterzones2(i)+secframes3));
                        elseif i == 1
                            bout2_3sec_120Same = cat(1,bout2_3sec_120Same,ztrial2(gcellind,enterzones2(i)-secframes3:enterzones2(i)+secframes3));
                        end
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes3
                            bout2_3sec_30Diff = cat(1,bout2_3sec_30Diff,ztrial2(gcellind,enterzones2(i)-secframes3:enterzones2(i)+secframes3));
                        elseif i == 1
                            bout2_3sec_30Diff = cat(1,bout2_3sec_30Diff,ztrial2(gcellind,enterzones2(i)-secframes3:enterzones2(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes3
                            bout2_3sec_120Diff = cat(1,bout2_3sec_120Diff,ztrial2(gcellind,enterzones2(i)-secframes3:enterzones2(i)+secframes3));
                        elseif i == 1
                            bout2_3sec_120Diff = cat(1,bout2_3sec_120Diff,ztrial2(gcellind,enterzones2(i)-secframes3:enterzones2(i)+secframes3));
                        end
                    end
                end
                % 2 pre 3 post interactions 
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes2
                            bout2_32sec_30Same = cat(1,bout2_32sec_30Same,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes3));
                        elseif i == 1
                            bout2_32sec_30Same = cat(1,bout2_32sec_30Same,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes2
                            bout2_32sec_120Same = cat(1,bout2_32sec_120Same,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes3));
                        elseif i == 1
                            bout2_32sec_120Same = cat(1,bout2_32sec_120Same,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes3));
                        end
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes2
                            bout2_32sec_30Diff = cat(1,bout2_32sec_30Diff,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes3));
                        elseif i == 1
                            bout2_32sec_30Diff = cat(1,bout2_32sec_30Diff,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzones2(i) - exitzones2(i-1) >= secframes2
                            bout2_32sec_120Diff = cat(1,bout2_32sec_120Diff,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes3));
                        elseif i == 1
                            bout2_32sec_120Diff = cat(1,bout2_32sec_120Diff,ztrial2(gcellind,enterzones2(i)-secframes2:enterzones2(i)+secframes3));
                        end
                    end
                end
            end
        end
    end
    % bout Hab1
    for i = 1: length(enterzonesHab1)
        if enterzonesHab1(i) > secframes2
            if exitzonesHab1(i) - enterzonesHab1(i) >= secframes2
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) >= secframes2
                            boutHab1_2sec_30Same = cat(1,boutHab1_2sec_30Same,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes2));
                        elseif i == 1
                            boutHab1_2sec_30Same = cat(1,boutHab1_2sec_30Same,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes2));
                        end
                        boutHab1_2sec_30Same_1h = boutcat_1h(enterzonesHab1,exitzonesHab1, secframes2,gcellind,zhab1,boutHab1_2sec_30Same_1h,i);
                        boutHab1_2sec_30Same_2h = boutcat_2h(enterzonesHab1,exitzonesHab1, secframes2,gcellind,zhab1,boutHab1_2sec_30Same_2h,i);
                    else
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) >= secframes2
                            boutHab1_2sec_120Same = cat(1,boutHab1_2sec_120Same,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes2));
                        elseif i == 1
                            boutHab1_2sec_120Same = cat(1,boutHab1_2sec_120Same,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes2));
                        end
                        boutHab1_2sec_120Same_1h = boutcat_1h(enterzonesHab1,exitzonesHab1, secframes2,gcellind,zhab1,boutHab1_2sec_120Same_1h,i);
                        boutHab1_2sec_120Same_2h = boutcat_2h(enterzonesHab1,exitzonesHab1, secframes2,gcellind,zhab1,boutHab1_2sec_120Same_2h,i);
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) >= secframes2
                            boutHab1_2sec_30Diff = cat(1,boutHab1_2sec_30Diff,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes2));
                        elseif i == 1
                            boutHab1_2sec_30Diff = cat(1,boutHab1_2sec_30Diff,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes2));
                        end
                        boutHab1_2sec_30Diff_1h = boutcat_1h(enterzonesHab1,exitzonesHab1, secframes2,gcellind,zhab1,boutHab1_2sec_30Diff_1h,i);
                        boutHab1_2sec_30Diff_2h = boutcat_2h(enterzonesHab1,exitzonesHab1, secframes2,gcellind,zhab1,boutHab1_2sec_30Diff_2h,i);
                    else
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) >= secframes2
                            boutHab1_2sec_120Diff = cat(1,boutHab1_2sec_120Diff,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes2));
                        elseif i == 1
                            boutHab1_2sec_120Diff = cat(1,boutHab1_2sec_120Diff,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes2));
                        end
                        boutHab1_2sec_120Diff_1h = boutcat_1h(enterzonesHab1,exitzonesHab1, secframes2,gcellind,zhab1,boutHab1_2sec_120Diff_1h,i);
                        boutHab1_2sec_120Diff_2h = boutcat_2h(enterzonesHab1,exitzonesHab1, secframes2,gcellind,zhab1,boutHab1_2sec_120Diff_2h,i);
                    end
                end
            end
        end
        if enterzonesHab1(i) > secframes3
            if exitzonesHab1(i) - enterzonesHab1(i) >= secframes3
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) > secframes3
                            boutHab1_3sec_30Same = cat(1,boutHab1_3sec_30Same,zhab1(gcellind,enterzonesHab1(i)-secframes3:enterzonesHab1(i)+secframes3));
                        elseif i == 1
                            boutHab1_3sec_30Same = cat(1,boutHab1_3sec_30Same,zhab1(gcellind,enterzonesHab1(i)-secframes3:enterzonesHab1(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) > secframes3
                            boutHab1_3sec_120Same = cat(1,boutHab1_3sec_120Same,zhab1(gcellind,enterzonesHab1(i)-secframes3:enterzonesHab1(i)+secframes3));
                        elseif i == 1
                            boutHab1_3sec_120Same = cat(1,boutHab1_3sec_120Same,zhab1(gcellind,enterzonesHab1(i)-secframes3:enterzonesHab1(i)+secframes3));
                        end
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) > secframes3
                            boutHab1_3sec_30Diff = cat(1,boutHab1_3sec_30Diff,zhab1(gcellind,enterzonesHab1(i)-secframes3:enterzonesHab1(i)+secframes3));
                        elseif i == 1
                            boutHab1_3sec_30Diff = cat(1,boutHab1_3sec_30Diff,zhab1(gcellind,enterzonesHab1(i)-secframes3:enterzonesHab1(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) > secframes3
                            boutHab1_3sec_120Diff = cat(1,boutHab1_3sec_120Diff,zhab1(gcellind,enterzonesHab1(i)-secframes3:enterzonesHab1(i)+secframes3));
                        elseif i == 1
                            boutHab1_3sec_120Diff = cat(1,boutHab1_3sec_120Diff,zhab1(gcellind,enterzonesHab1(i)-secframes3:enterzonesHab1(i)+secframes3));
                        end
                    end
                end
                % 2 pre 3 post interactions 
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) > secframes2
                            boutHab1_32sec_30Same = cat(1,boutHab1_32sec_30Same,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes3));
                        elseif i == 1
                            boutHab1_32sec_30Same = cat(1,boutHab1_32sec_30Same,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) > secframes2
                            boutHab1_32sec_120Same = cat(1,boutHab1_32sec_120Same,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes3));
                        elseif i == 1
                            boutHab1_32sec_120Same = cat(1,boutHab1_32sec_120Same,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes3));
                        end
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) > secframes2
                            boutHab1_32sec_30Diff = cat(1,boutHab1_32sec_30Diff,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes3));
                        elseif i == 1
                            boutHab1_32sec_30Diff = cat(1,boutHab1_32sec_30Diff,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzonesHab1(i) - exitzonesHab1(i-1) > secframes2
                            boutHab1_32sec_120Diff = cat(1,boutHab1_32sec_120Diff,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes3));
                        elseif i == 1
                            boutHab1_32sec_120Diff = cat(1,boutHab1_32sec_120Diff,zhab1(gcellind,enterzonesHab1(i)-secframes2:enterzonesHab1(i)+secframes3));
                        end
                    end
                end
            end            
        end
    end
    % bout Hab2
    for i = 1: length(enterzonesHab2)
        if enterzonesHab2(i) > secframes2
            if exitzonesHab2(i) - enterzonesHab2(i) >= secframes2
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) >= secframes2
                            boutHab2_2sec_30Same = cat(1,boutHab2_2sec_30Same,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes2));
                        elseif i == 1
                            boutHab2_2sec_30Same = cat(1,boutHab2_2sec_30Same,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes2));
                        end
                        boutHab2_2sec_30Same_1h = boutcat_1h(enterzonesHab2,exitzonesHab2, secframes2,gcellind,zhab2,boutHab2_2sec_30Same_1h,i);
                        boutHab2_2sec_30Same_2h = boutcat_2h(enterzonesHab2,exitzonesHab2, secframes2,gcellind,zhab2,boutHab2_2sec_30Same_2h,i);
                    else
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) >= secframes2
                            boutHab2_2sec_120Same = cat(1,boutHab2_2sec_120Same,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes2));
                        elseif i == 1
                            boutHab2_2sec_120Same = cat(1,boutHab2_2sec_120Same,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes2));
                        end
                        boutHab2_2sec_120Same_1h = boutcat_1h(enterzonesHab2,exitzonesHab2, secframes2,gcellind,zhab2,boutHab2_2sec_120Same_1h,i);
                        boutHab2_2sec_120Same_2h = boutcat_2h(enterzonesHab2,exitzonesHab2, secframes2,gcellind,zhab2,boutHab2_2sec_120Same_2h,i);
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) >= secframes2
                            boutHab2_2sec_30Diff = cat(1,boutHab2_2sec_30Diff,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes2));
                        elseif i == 1
                            boutHab2_2sec_30Diff = cat(1,boutHab2_2sec_30Diff,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes2));
                        end
                        boutHab2_2sec_30Diff_1h = boutcat_1h(enterzonesHab2,exitzonesHab2, secframes2,gcellind,zhab2,boutHab2_2sec_30Diff_1h,i);
                        boutHab2_2sec_30Diff_2h = boutcat_2h(enterzonesHab2,exitzonesHab2, secframes2,gcellind,zhab2,boutHab2_2sec_30Diff_2h,i);
                    else
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) >= secframes2
                            boutHab2_2sec_120Diff = cat(1,boutHab2_2sec_120Diff,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes2));
                        elseif i == 1
                            boutHab2_2sec_120Diff = cat(1,boutHab2_2sec_120Diff,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes2));
                        end
                        boutHab2_2sec_120Diff_1h = boutcat_1h(enterzonesHab2,exitzonesHab2, secframes2,gcellind,zhab2,boutHab2_2sec_120Diff_1h,i);
                        boutHab2_2sec_120Diff_2h = boutcat_2h(enterzonesHab2,exitzonesHab2, secframes2,gcellind,zhab2,boutHab2_2sec_120Diff_2h,i);
                    end
                end
            end
        end
        if enterzonesHab2(i) > secframes3
            if exitzonesHab2(i) - enterzonesHab2(i) >= secframes3
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) > secframes3
                            boutHab2_3sec_30Same = cat(1,boutHab2_3sec_30Same,zhab2(gcellind,enterzonesHab2(i)-secframes3:enterzonesHab2(i)+secframes3));
                        elseif i == 1
                            boutHab2_3sec_30Same = cat(1,boutHab2_3sec_30Same,zhab2(gcellind,enterzonesHab2(i)-secframes3:enterzonesHab2(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) > secframes3
                            boutHab2_3sec_120Same = cat(1,boutHab2_3sec_120Same,zhab2(gcellind,enterzonesHab2(i)-secframes3:enterzonesHab2(i)+secframes3));
                        elseif i == 1
                            boutHab2_3sec_120Same = cat(1,boutHab2_3sec_120Same,zhab2(gcellind,enterzonesHab2(i)-secframes3:enterzonesHab2(i)+secframes3));
                        end
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) > secframes3
                            boutHab2_3sec_30Diff = cat(1,boutHab2_3sec_30Diff,zhab2(gcellind,enterzonesHab2(i)-secframes3:enterzonesHab2(i)+secframes3));
                        elseif i == 1
                            boutHab2_3sec_30Diff = cat(1,boutHab2_3sec_30Diff,zhab2(gcellind,enterzonesHab2(i)-secframes3:enterzonesHab2(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) > secframes3
                            boutHab2_3sec_120Diff = cat(1,boutHab2_3sec_120Diff,zhab2(gcellind,enterzonesHab2(i)-secframes3:enterzonesHab2(i)+secframes3));
                        elseif i == 1
                            boutHab2_3sec_120Diff = cat(1,boutHab2_3sec_120Diff,zhab2(gcellind,enterzonesHab2(i)-secframes3:enterzonesHab2(i)+secframes3));
                        end
                    end
                end
            end
            %2 sec pre and 3 sec post interaction
            if exitzonesHab2(i) - enterzonesHab2(i) >= secframes3
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) > secframes2
                            boutHab2_32sec_30Same = cat(1,boutHab2_32sec_30Same,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes3));
                        elseif i == 1
                            boutHab2_32sec_30Same = cat(1,boutHab2_32sec_30Same,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) > secframes2
                            boutHab2_32sec_120Same = cat(1,boutHab2_32sec_120Same,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes3));
                        elseif i == 1
                            boutHab2_32sec_120Same = cat(1,boutHab2_32sec_120Same,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes3));
                        end
                    end
                else
                    if SITnovelty_mins{SITseshnum,2} == 30
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) > secframes2
                            boutHab2_32sec_30Diff = cat(1,boutHab2_32sec_30Diff,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes3));
                        elseif i == 1
                            boutHab2_32sec_30Diff = cat(1,boutHab2_32sec_30Diff,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes3));
                        end
                    else
                        if i > 1 && enterzonesHab2(i) - exitzonesHab2(i-1) > secframes2
                            boutHab2_32sec_120Diff = cat(1,boutHab2_32sec_120Diff,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes3));
                        elseif i == 1
                            boutHab2_32sec_120Diff = cat(1,boutHab2_32sec_120Diff,zhab2(gcellind,enterzonesHab2(i)-secframes2:enterzonesHab2(i)+secframes3));
                        end
                    end
                end
            end
        end
    end
    
    padval1 = 50 - mod(length(catActivity1(1,:)),timebin2);
    padval2 = 50 - mod(length(catActivity2(1,:)),timebin2);
    padvalHab1 = 50 - mod(length(catActivityHab1(1,:)),timebin2);
    padvalHab2 = 50 - mod(length(catActivityHab2(1,:)),timebin2);
    
    catActivity1pad = padarray(catActivity1,[0,padval1],nan,'post');
    catActivity2pad = padarray(catActivity2,[0,padval2],nan,'post');
    catActivityHab1pad = padarray(catActivityHab1,[0,padvalHab1],nan,'post');
    catActivityHab2pad = padarray(catActivityHab2,[0,padvalHab2],nan,'post');
    
    scalingfact1 = round(length(catActivity1pad(1,:))/50);
    scalingfact2 = round(length(catActivity2pad(1,:))/50);
    scalingfacthab1 = round(length(catActivityHab1pad(1,:))/50);
    scalingfacthab2 = round(length(catActivityHab2pad(1,:))/50);
    
    blockSize = [1, 2];
    meanFilterFunction = @(theBlockStructure) nansum(theBlockStructure.data(:));
%     blockAveragedDownSignal = blockproc(PulseRateF, blockSize, meanFilterFunction);
    
    for i = 1 : celltotal
        ChronoActivity1(i,:) = blockproc(catActivity1pad(i,:),[1,scalingfact1],meanFilterFunction)./scalingfact1;
        ChronoActivity2(i,:) = blockproc(catActivity2pad(i,:),[1,scalingfact2],meanFilterFunction)./scalingfact2;
        ChronoActivityHab1(i,:) = blockproc(catActivityHab1pad(i,:),[1,scalingfacthab1],meanFilterFunction)./scalingfacthab1;
        ChronoActivityHab2(i,:) = blockproc(catActivityHab2pad(i,:),[1,scalingfacthab2],meanFilterFunction)./scalingfacthab2;
    end
    
    chronoall{SITseshnum,1} = ChronoActivityHab1;
    chronoall{SITseshnum,2} = ChronoActivity1;
    chronoall{SITseshnum,3} = ChronoActivityHab2;
    chronoall{SITseshnum,4} = ChronoActivity2;
    
    if contains(SITnovelty_mins{SITseshnum,1},'Same')
        if SITnovelty_mins{SITseshnum,2} == 30 
            if isempty(ChronoAllcells30Hab1Same)
                ChronoAllcells30Hab1Same = ChronoActivityHab1;
                ChronoAllcells30Trial1Same = ChronoActivity1;
                ChronoAllcells30Hab2Same = ChronoActivityHab2;
                ChronoAllcells30Trial2Same = ChronoActivity2;
            else
                ChronoAllcells30Hab1Same = cat(1,ChronoAllcells30Hab1Same,ChronoActivityHab1);
                ChronoAllcells30Trial1Same = cat(1,ChronoAllcells30Trial1Same,ChronoActivity1);
                ChronoAllcells30Hab2Same = cat(1,ChronoAllcells30Hab2Same,ChronoActivityHab2);
                ChronoAllcells30Trial2Same = cat(1,ChronoAllcells30Trial2Same,ChronoActivity2);
            end
        else
            if isempty(ChronoAllcells120Hab1Same)
                ChronoAllcells120Hab1Same = ChronoActivityHab1;
                ChronoAllcells120Trial1Same = ChronoActivity1;
                ChronoAllcells120Hab2Same = ChronoActivityHab2;
                ChronoAllcells120Trial2Same = ChronoActivity2;
            else
                ChronoAllcells120Hab1Same = cat(1,ChronoAllcells120Hab1Same,ChronoActivityHab1);
                ChronoAllcells120Trial1Same = cat(1,ChronoAllcells120Trial1Same,ChronoActivity1);
                ChronoAllcells120Hab2Same = cat(1,ChronoAllcells120Hab2Same,ChronoActivityHab2);
                ChronoAllcells120Trial2Same = cat(1,ChronoAllcells120Trial2Same,ChronoActivity2);
            end
        end
    else
        if SITnovelty_mins{SITseshnum,2} == 30
            if isempty(ChronoAllcells30Hab1Diff)
                ChronoAllcells30Hab1Diff = ChronoActivityHab1;
                ChronoAllcells30Trial1Diff = ChronoActivity1;
                ChronoAllcells30Hab2Diff = ChronoActivityHab2;
                ChronoAllcells30Trial2Diff = ChronoActivity2;
            else
                ChronoAllcells30Hab1Diff = cat(1,ChronoAllcells30Hab1Diff,ChronoActivityHab1);
                ChronoAllcells30Trial1Diff = cat(1,ChronoAllcells30Trial1Diff,ChronoActivity1);
                ChronoAllcells30Hab2Diff = cat(1,ChronoAllcells30Hab2Diff,ChronoActivityHab2);
                ChronoAllcells30Trial2Diff = cat(1,ChronoAllcells30Trial2Diff,ChronoActivity2);
            end
        else
            if isempty(ChronoAllcells120Hab1Diff)
                ChronoAllcells120Hab1Diff = ChronoActivityHab1;
                ChronoAllcells120Trial1Diff = ChronoActivity1;
                ChronoAllcells120Hab2Diff = ChronoActivityHab2;
                ChronoAllcells120Trial2Diff = ChronoActivity2;
            else
                ChronoAllcells120Hab1Diff = cat(1,ChronoAllcells120Hab1Diff,ChronoActivityHab1);
                ChronoAllcells120Trial1Diff = cat(1,ChronoAllcells120Trial1Diff,ChronoActivity1);
                ChronoAllcells120Hab2Diff = cat(1,ChronoAllcells120Hab2Diff,ChronoActivityHab2);
                ChronoAllcells120Trial2Diff = cat(1,ChronoAllcells120Trial2Diff,ChronoActivity2);
            end
        end
    end
    
    outztrial1 = ztrial1;
    outztrial2 = ztrial2;
    outzhab1 = zhab1;
    outzhab2 = zhab2;
    
    for i = length(enterzones1) : -1 : 1
        outztrial1(:,enterzones1(i) : exitzones1(i)) = [];
    end
    for i = length(enterzones2) : -1 : 1
        outztrial2(:,enterzones2(i) : exitzones2(i)) = [];
    end
    for i = length(enterzonesHab1) : -1 : 1
        outzhab1(:,enterzonesHab1(i) : exitzonesHab1(i)) = [];
    end
    for i = length(enterzonesHab2) : -1 : 1
        outzhab2(:,enterzonesHab2(i) : exitzonesHab2(i)) = [];
    end
    
    outZone1 = sum(outztrial1(gcellind,:),2)./length(outztrial1(1,:));
    outZone2 = sum(outztrial2(gcellind,:),2)./length(outztrial2(1,:));
    outZoneHab1 = sum(outzhab1(gcellind,:),2)./length(outzhab1(1,:));
    outZoneHab2 = sum(outzhab2(gcellind,:),2)./length(outzhab2(1,:));
    
    shuffleitterations = 100;
    stdInteractionActivityTrial1 = std(zoneActivityRaw1,0,2);
    stdInteractionActivityTrial2 = std(zoneActivityRaw2,0,2);
    stdInteractionActivityHab1 = std(zoneActivityHabRaw1,0,2);
    stdInteractionActivityHab2 = std(zoneActivityHabRaw2,0,2);
    shuffledInteractions1 = ShuffleSITinterations(trial1,enterzones1,exitzones1,shuffleitterations);
    shuffledInteractions2 = ShuffleSITinterations(trial2,enterzones2,exitzones2,shuffleitterations);
    shuffledInteractionsHab1 = ShuffleSITinterations(hab1,enterzonesHab1,exitzonesHab1,shuffleitterations);
    shuffledInteractionsHab2 = ShuffleSITinterations(hab2,enterzonesHab2,exitzonesHab2,shuffleitterations);
    
    shuffledInteractionsZ1 = ShuffleSITinterations(ztrials,enterzones1+length(zhab1(1,:)),exitzones1+length(zhab1(1,:)),shuffleitterations);
    shuffledInteractionsZ2 = ShuffleSITinterations(ztrials,enterzones2+length(ztrial1(1,:))+length(zhab1(1,:))+length(zhab2(1,:)),exitzones2+length(ztrial1(1,:))+length(zhab1(1,:))+length(zhab2(1,:)),shuffleitterations);
    shuffledInteractionsZHab1 = ShuffleSITinterations(ztrials,enterzonesHab1,exitzonesHab1,shuffleitterations);
    shuffledInteractionsZHab2 = ShuffleSITinterations(ztrials,enterzonesHab2+length(ztrial1(1,:))+length(zhab1(1,:)),exitzonesHab2+length(ztrial1(1,:))+length(zhab1(1,:)),shuffleitterations);
    
    
    shuffledInteractionsZ1_2 = ShuffleSITinterations(ztrials_2,enterzones1+length(zhab1(1,4501:end)),exitzones1+length(zhab1(1,4501:end)),shuffleitterations);
    shuffledInteractionsZ2_2 = ShuffleSITinterations(ztrials_2,enterzones2+length(ztrial1(1,:))+length(zhab1(1,4501:end))+length(zhab2(1,4501:end)),exitzones2+length(ztrial1(1,:))+length(zhab1(1,4501:end))+length(zhab2(1,4501:end)),shuffleitterations);
    shuffledInteractionsZHab1_2 = ShuffleSITinterations(ztrials_2,enterzonesHab1_2,exitzonesHab1_2,shuffleitterations);
    shuffledInteractionsZHab2_2 = ShuffleSITinterations(ztrials_2,enterzonesHab2_2+length(ztrial1(1,:))+length(zhab1(1,4501:end)),exitzonesHab2_2+length(ztrial1(1,:))+length(zhab1(1,4501:end)),shuffleitterations);
    
    save('InteractionCalcium','catActivity1','catActivity1Raw','catActivity2','catActivity2Raw','catActivityHab1','catActivityHab1Raw','catActivityHab2','catActivityHab2Raw')
    save('MeanFiringRateInteractions.mat','zoneActivityRaw1','zoneActivityRaw2','zoneActivityHabRaw1','zoneActivityHabRaw2','stdInteractionActivityTrial1','stdInteractionActivityTrial2','stdInteractionActivityHab1','stdInteractionActivityHab2',...
        'zoneActivity1','zoneActivity2','zoneActivityHab1','zoneActivityHab2',...
        'zoneActivityHab1_2','zoneActivityHab2_2')
    save('ShuffledMeanFiringRateInteractions.mat','shuffledInteractions1','shuffledInteractions2','shuffledInteractionsHab1','shuffledInteractionsHab2',...
        'shuffledInteractionsZ1','shuffledInteractionsZ2','shuffledInteractionsZHab1','shuffledInteractionsZHab2',...
        'shuffledInteractionsZ1_2','shuffledInteractionsZ2_2','shuffledInteractionsZHab1_2','shuffledInteractionsZHab2_2')
    
    save('Figure4b50binData','catActivity1','catActivity2','catActivityHab1','catActivityHab2','catActivity1Raw','catActivity2Raw','catActivityHab1Raw','catActivityHab2Raw','ChronoActivity1','ChronoActivity2','ChronoActivityHab1','ChronoActivityHab2')
    
    meanZoneActivity1 = zeros(celltotal,timebinnum);
    meanZoneActivity2 = zeros(celltotal,timebinnum);
    meanZoneActivityHab1 = zeros(celltotal,timebinnum);
    meanZoneActivityHab2 = zeros(celltotal,timebinnum);
    for i = 1 : timebinnum
        for j = 1 : celltotal
            meanZoneActivity1(j,i) = mean(zoneActivitysplit1(j,:,i));
            meanZoneActivity2(j,i) = mean(zoneActivitysplit2(j,:,i));
            meanZoneActivityHab1(j,i) = mean(zoneActivitysplitHab1(j,:,i));
            meanZoneActivityHab2(j,i) = mean(zoneActivitysplitHab2(j,:,i));
        end
    end
    meanacttot1 = mean(meanZoneActivity1,2);
    meanacttot2 = mean(meanZoneActivity2,2);
    meanacttotHab1 = mean(meanZoneActivityHab1,2);
    meanacttotHab2 = mean(meanZoneActivityHab2,2);
    
    maxvals = cat(1,max(max(meanZoneActivity1)), max(max(meanZoneActivity2)),max(max(meanZoneActivityHab1)), max(max(meanZoneActivityHab2)));
    maxind = find(maxvals == max(maxvals));
    maxmean = cat(1,max(max(meanacttot1)), max(max(meanacttot2)),max(max(meanacttotHab1)), max(max(meanacttotHab2)));
    maxmeanind = find(maxmean == max(maxmean));
    
    if maxind == 1
        maxZone = max(max(meanZoneActivity1));
    elseif maxind == 2
        maxZone = max(max(meanZoneActivity2));
    elseif maxind == 3
        maxZone = max(max(meanZoneActivityHab1));
    elseif maxind == 4
        maxZone = max(max(meanZoneActivityHab2));
    end
    
    if maxmeanind == 1
        maxMean = max(max(meanacttot1));
    elseif maxmeanind == 2
        maxMean = max(max(meanacttot2));
    elseif maxmeanind == 3
        maxMean = max(max(meanacttotHab1));
    elseif maxmeanind == 4
        maxMean = max(max(meanacttotHab2));
    end
    
    maxchrono = max(max(cat(1,ChronoActivity1,ChronoActivity2,ChronoActivityHab1,ChronoActivityHab2)));
    
    [~, indActivity] = sort(meanacttot1,'descend');
    
    figure
    subplot(1,4,1)
    imagesc(meanZoneActivityHab1(indActivity,:))
    caxis([0 maxZone])
    colorbar
    title('Hab1')
    subplot(1,4,2)
    imagesc(meanZoneActivity1(indActivity,:))
    caxis([0 maxZone])
    colorbar
    title('Trial1')
    subplot(1,4,3)
    imagesc(meanZoneActivityHab2(indActivity,:))
    caxis([0 maxZone])
    title('Hab2')
    colorbar
    subplot(1,4,4)
    imagesc(meanZoneActivity2(indActivity,:))
    caxis([0 maxZone])
    title('Trial2')
    colorbar
    sgtitle([seshdate ' ' mousename ' ' SITnovelty_mins{SITseshnum,1} ' ' num2str(SITnovelty_mins{SITseshnum,2}) 'mins'])
    
    saveas(gcf,['Figure4b'])
    saveas(gcf,['Figure4b.jpeg'])
    
    figure
    subplot(1,4,1)
    imagesc(meanacttotHab1(indActivity,:))
    caxis([0 maxMean])
    colorbar
    title('Hab1')
    subplot(1,4,2)
    imagesc(meanacttot1(indActivity,:))
    caxis([0 maxMean])
    colorbar
    title('Trial1')
    subplot(1,4,3)
    imagesc(meanacttotHab2(indActivity,:))
    caxis([0 maxMean])
    title('Hab2')
    colorbar
    subplot(1,4,4)
    imagesc(meanacttot2(indActivity,:))
    caxis([0 maxMean])
    title('Trial2')
    colorbar
    sgtitle([seshdate ' ' mousename ' ' SITnovelty_mins{SITseshnum,1} ' ' num2str(SITnovelty_mins{SITseshnum,2}) 'mins'])
    
    saveas(gcf,['MeanInteractionFiringRate'])
    saveas(gcf,['MeanInteractionFiringRate.jpeg'])
    save('FiringRateSortedIndex.mat','indActivity','gcellind')
    
    figure
    subplot(1,8,1)
    imagesc(meanacttotHab1(indActivity,:))
    caxis([0 maxMean])
    colorbar
    title('In Hab1')
    subplot(1,8,2)
    imagesc(outZoneHab1(indActivity,:))
    caxis([0 maxMean])
    colorbar
    title('Out Hab1')
    subplot(1,8,3)
    imagesc(meanacttot1(indActivity,:))
    caxis([0 maxMean])
    colorbar
    title('In Trial1')
    subplot(1,8,4)
    imagesc(outZone1(indActivity,:))
    caxis([0 maxMean])
    colorbar
    title('Out Trial1')
    subplot(1,8,5)
    imagesc(meanacttotHab2(indActivity,:))
    caxis([0 maxMean])
    title('In Hab2')
    colorbar
    subplot(1,8,6)
    imagesc(outZoneHab2(indActivity,:))
    caxis([0 maxMean])
    title('Out Hab2')
    colorbar
    subplot(1,8,7)
    imagesc(meanacttot2(indActivity,:))
    caxis([0 maxMean])
    title('In Trial2')
    colorbar
    subplot(1,8,8)
    imagesc(outZone2(indActivity,:))
    caxis([0 maxMean])
    title('Out Trial2')
    colorbar
    sgtitle([seshdate ' ' mousename ' ' SITnovelty_mins{SITseshnum,1} ' ' num2str(SITnovelty_mins{SITseshnum,2}) 'mins'])
    
    saveas(gcf,['MeanInteraction&OutInteractionFiringRate'])
    saveas(gcf,['MeanInteraction&OutInteractionFiringRate.jpeg'])
    save('CellIndInfo','celltotal','gcellind','indActivity')
    save('outinteraction.mat','outZone1','outZone2','outZoneHab1','outZoneHab2')
    save('zdataUncut.mat','ztrial1','ztrial2','zhab1','zhab2')
    
    
    figure
    subplot(1,4,1)
    imagesc(ChronoActivityHab1(indActivity,:))
    caxis([0 maxchrono])
    colorbar
    ylabel('Cell Number')
    title('Hab1')
    subplot(1,4,2)
    imagesc(ChronoActivity1(indActivity,:))
    caxis([0 maxchrono])
    colorbar
    title('Trial1')
    subplot(1,4,3)
    imagesc(ChronoActivityHab2(indActivity,:))
    caxis([0 maxchrono])
    colorbar
    title('Hab2')
    subplot(1,4,4)
    imagesc(ChronoActivity2(indActivity,:))
    caxis([0 maxchrono])
    colorbar
    title('Trial2')
    sgtitle([seshdate ' ' mousename ' ' SITnovelty_mins{SITseshnum,1} ' ' num2str(SITnovelty_mins{SITseshnum,2}) 'mins'])
    saveas(gcf,['Figure4b_50bins'])
    saveas(gcf,['Figure4b_50bins.jpeg'])
    
    %}
    
    if ~isempty(find(find(strcmp(SITproximityData.micenames,mousename)) < SITseshnum))
        previousSesh = find(strcmp(SITproximityData.micenames,mousename));
        
        folderparts = strsplit(SITproximityData.folderpaths{1,previousSesh(1)},'\');
        datepartsPrev = strsplit(folderparts{end-3},' ');
        seshdatePrev = datepartsPrev{2};
        previnteractions = load([cdold '\' seshdatePrev '\' mousename '\MeanFiringRateInteractions.mat']);
        prevshuffling = load([cdold '\' seshdatePrev '\' mousename '\ShuffledMeanFiringRateInteractions.mat']);
        prevCellinds = load([cdold '\' seshdatePrev '\' mousename '\FiringRateSortedIndex.mat']);
        prev4b = load([cdold '\' seshdatePrev '\' mousename '\Figure4b50binData.mat']);
        
        shuffledCutoff1 = permute(prctile(shuffledInteractions1(gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoff2 = permute(prctile(shuffledInteractions2(gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffHab1 = permute(prctile(shuffledInteractionsHab1(gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffHab2 = permute(prctile(shuffledInteractionsHab2(gcellind,:,:),95,2),[1 3 2]);
        
        shuffledCutoff1prev = permute(prctile(prevshuffling.shuffledInteractions1(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoff2prev = permute(prctile(prevshuffling.shuffledInteractions2(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffHab1prev = permute(prctile(prevshuffling.shuffledInteractionsHab1(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffHab2prev = permute(prctile(prevshuffling.shuffledInteractionsHab2(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        
        PassedShuffledCutoff1 = zoneActivityRaw1 >= shuffledCutoff1;
        PassedShuffledCutoff2 = zoneActivityRaw2 >= shuffledCutoff2;
        PassedShuffledCutoffHab1 = zoneActivityHabRaw1 >= shuffledCutoffHab1;
        PassedShuffledCutoffHab2 = zoneActivityHabRaw2 >= shuffledCutoffHab2;
        
        PassedShuffledCutoff1prev = previnteractions.zoneActivity1 > shuffledCutoff1prev;
        PassedShuffledCutoff2prev = previnteractions.zoneActivity2 > shuffledCutoff2prev;
        PassedShuffledCutoffHab1prev = previnteractions.zoneActivityHab1 > shuffledCutoffHab1prev;
        PassedShuffledCutoffHab2prev = previnteractions.zoneActivityHab2 > shuffledCutoffHab2prev;
        
        p1 = find(sum(PassedShuffledCutoff1,2)/length(PassedShuffledCutoff1(1,:))>=0.5);
        p2 = find(sum(PassedShuffledCutoff2,2)/length(PassedShuffledCutoff2(1,:))>=0.5);
        ph1 = find(sum(PassedShuffledCutoffHab1,2)/length(PassedShuffledCutoffHab1(1,:))>=0.5);
        ph2 = find(sum(PassedShuffledCutoffHab2,2)/length(PassedShuffledCutoffHab2(1,:))>=0.5);
        
        p1prev = find(sum(PassedShuffledCutoff1prev,2)/length(PassedShuffledCutoff1prev(1,:))>=0.5);
        p2prev = find(sum(PassedShuffledCutoff2prev,2)/length(PassedShuffledCutoff2prev(1,:))>=0.5);
        ph1prev = find(sum(PassedShuffledCutoffHab1prev,2)/length(PassedShuffledCutoffHab1prev(1,:))>=0.5);
        ph2prev = find(sum(PassedShuffledCutoffHab2prev,2)/length(PassedShuffledCutoffHab2prev(1,:))>=0.5);
                        
        shuffledCutoffZ1 = permute(prctile(shuffledInteractionsZ1(gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZ2 = permute(prctile(shuffledInteractionsZ2(gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZHab1 = permute(prctile(shuffledInteractionsZHab1(gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZHab2 = permute(prctile(shuffledInteractionsZHab2(gcellind,:,:),95,2),[1 3 2]);
        
        shuffledCutoffZ1prev = permute(prctile(prevshuffling.shuffledInteractionsZ1(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZ2prev = permute(prctile(prevshuffling.shuffledInteractionsZ2(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZHab1prev = permute(prctile(prevshuffling.shuffledInteractionsZHab1(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZHab2prev = permute(prctile(prevshuffling.shuffledInteractionsZHab2(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        
        PassedShuffledCutoffZ1 = zoneActivity1 >= shuffledCutoffZ1;
        PassedShuffledCutoffZ2 = zoneActivity2 >= shuffledCutoffZ2;
        PassedShuffledCutoffZHab1 = zoneActivityHab1 >= shuffledCutoffZHab1;
        PassedShuffledCutoffZHab2 = zoneActivityHab2 >= shuffledCutoffZHab2;
        
        PassedShuffledCutoffZ1prev = previnteractions.zoneActivity1 >= shuffledCutoffZ1prev;
        PassedShuffledCutoffZ2prev = previnteractions.zoneActivity2 >= shuffledCutoffZ2prev;
        PassedShuffledCutoffZHab1prev = previnteractions.zoneActivityHab1 >= shuffledCutoffZHab1prev;
        PassedShuffledCutoffZHab2prev = previnteractions.zoneActivityHab2 >= shuffledCutoffZHab2prev;  
        
        p1z = find(sum(PassedShuffledCutoffZ1,2));
        p2z = find(sum(PassedShuffledCutoffZ2,2));
        ph1z = find(sum(PassedShuffledCutoffZHab1,2));
        ph2z = find(sum(PassedShuffledCutoffZHab2,2));
        
        p1z_not = find(sum(PassedShuffledCutoffZ1,2)== 0);
        p2z_not = find(sum(PassedShuffledCutoffZ2,2)== 0);
        ph1z_not = find(sum(PassedShuffledCutoffZHab1,2)== 0);
        ph2z_not = find(sum(PassedShuffledCutoffZHab2,2)== 0);
        
        p1prevz = find(sum(PassedShuffledCutoffZ1prev,2));
        p2prevz = find(sum(PassedShuffledCutoffZ2prev,2));
        ph1prevz = find(sum(PassedShuffledCutoffZHab1prev,2));
        ph2prevz = find(sum(PassedShuffledCutoffZHab2prev,2));
        
        p1prevz_not = find(sum(PassedShuffledCutoffZ1prev,2)==0);
        p2prevz_not = find(sum(PassedShuffledCutoffZ2prev,2)==0);
        ph1prevz_not = find(sum(PassedShuffledCutoffZHab1prev,2)==0);
        ph2prevz_not = find(sum(PassedShuffledCutoffZHab2prev,2)==0);
        
        shuffledCutoffZ1_2 = permute(prctile(shuffledInteractionsZ1_2(gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZ2_2 = permute(prctile(shuffledInteractionsZ2_2(gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZHab1_2 = permute(prctile(shuffledInteractionsZHab1_2(gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZHab2_2 = permute(prctile(shuffledInteractionsZHab2_2(gcellind,:,:),95,2),[1 3 2]);
        
        shuffledCutoffZ1prev_2 = permute(prctile(prevshuffling.shuffledInteractionsZ1_2(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZ2prev_2 = permute(prctile(prevshuffling.shuffledInteractionsZ2_2(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZHab1prev_2 = permute(prctile(prevshuffling.shuffledInteractionsZHab1_2(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        shuffledCutoffZHab2prev_2 = permute(prctile(prevshuffling.shuffledInteractionsZHab2_2(prevCellinds.gcellind,:,:),95,2),[1 3 2]);
        
        PassedShuffledCutoffZ1_2 = zoneActivity1 > shuffledCutoffZ1_2;
        PassedShuffledCutoffZ2_2 = zoneActivity2 > shuffledCutoffZ2_2;
        PassedShuffledCutoffZHab1_2 = zoneActivityHab1_2 > shuffledCutoffZHab1_2;
        PassedShuffledCutoffZHab2_2 = zoneActivityHab2_2 > shuffledCutoffZHab2_2;
        
        PassedShuffledCutoffZ1prev_2 = previnteractions.zoneActivity1 > shuffledCutoffZ1prev_2;
        PassedShuffledCutoffZ2prev_2 = previnteractions.zoneActivity2 > shuffledCutoffZ2prev_2;
        PassedShuffledCutoffZHab1prev_2 = previnteractions.zoneActivityHab1_2 > shuffledCutoffZHab1prev_2;
        PassedShuffledCutoffZHab2prev_2 = previnteractions.zoneActivityHab2_2 > shuffledCutoffZHab2prev_2;
        
        p1z_CR = find(sum(PassedShuffledCutoffZ1,2)/length(PassedShuffledCutoffZ1(1,:)));
        p2z_CR = find(sum(PassedShuffledCutoffZ2,2)/length(PassedShuffledCutoffZ2(1,:)));
        ph1z_CR = find(sum(PassedShuffledCutoffZHab1,2)/length(PassedShuffledCutoffZHab1(1,:)));
        ph2z_CR = find(sum(PassedShuffledCutoffZHab2,2)/length(PassedShuffledCutoffZHab2(1,:)));
        
        p1prevz_CR = find(sum(PassedShuffledCutoffZ1prev,2))/length(PassedShuffledCutoffZ1prev(1,:));
        p2prevz_CR = find(sum(PassedShuffledCutoffZ2prev,2))/length(PassedShuffledCutoffZ2prev(1,:));
        ph1prevz_CR = find(sum(PassedShuffledCutoffZHab1prev,2))/length(PassedShuffledCutoffZHab1prev(1,:));
        ph2prevz_CR = find(sum(PassedShuffledCutoffZHab2prev,2))/length(PassedShuffledCutoffZHab2prev(1,:));                
        
        figure
        if contains(SITnovelty_mins{SITseshnum,1},'Different')
            subplot(1,2,2)
            plot([1,2],[mean(zoneActivityRaw1,2),mean(zoneActivityRaw2,2)])
            xlim([0.8 2.2])
            title([seshdate ' ' SITnovelty_mins{SITseshnum,1} ' Mouse'])
            xticks([1,2])
            xticklabels({'Trial1', 'Trial2'})
            subplot(1,2,1)
            plot([1,2],[mean(previnteractions.zoneActivityRaw1,2),mean(previnteractions.zoneActivityRaw2,2)])
            xlim([0.8 2.2])
            title([seshdatePrev ' ' SITnovelty_mins{previousSesh(1),1} ' Mouse'])
            xticks([1,2])
            xticklabels({'Trial1', 'Trial2'})
        else
            subplot(1,2,1)
            plot([1,2],[mean(zoneActivityRaw1,2),mean(zoneActivityRaw2,2)])
            xlim([0.8 2.2])
            title([seshdate ' ' SITnovelty_mins{SITseshnum,1} ' Mouse'])
            xticks([1,2])
            xticklabels({'Trial1', 'Trial2'})
            subplot(1,2,2)
            plot([1,2],[mean(previnteractions.zoneActivityRaw1,2),mean(previnteractions.zoneActivityRaw2,2)])
            xlim([0.8 2.2])
            title([seshdatePrev ' ' SITnovelty_mins{previousSesh(1),1} ' Mouse'])
            xticks([1,2])
            xticklabels({'Trial1', 'Trial2'})
        end
        sgtitle([mousename ' Cell level Interaction Mean Firing Rate across days' num2str(SITnovelty_mins{SITseshnum,2}) 'mins'])
        saveas(gcf,['CellInteractionMeanFiringAcrossTrials'])
        saveas(gcf,['CellInteractionMeanFiringAcrossTrials.jpeg'])
        
        figure
        Interactions1mean = mean(mean(zoneActivityRaw1,2));
        stdmean1 = std(mean(zoneActivityRaw1,2))/sqrt(length(zoneActivityRaw1));
        Interactions2mean = mean(mean(zoneActivityRaw2,2));
        stdmean2 = std(mean(zoneActivityRaw2,2))/sqrt(length(zoneActivityRaw2));
        Interactions1Prevmean = mean(mean(previnteractions.zoneActivityRaw1,2));
        stdPrevmean1 = std(mean(previnteractions.zoneActivityRaw1,2))/sqrt(length(previnteractions.zoneActivityRaw1));
        Interactions2Prevmean = mean(mean(previnteractions.zoneActivityRaw2,2));
        stdPrevmean2 = std(mean(previnteractions.zoneActivityRaw2,2))/sqrt(length(previnteractions.zoneActivityRaw2));
        
        if contains(SITnovelty_mins{SITseshnum,1},'Same')
            scatter([1;2;3;4],[Interactions1mean;Interactions2mean;Interactions1Prevmean;Interactions2Prevmean])
            errorbar([1,2,3,4],[Interactions1mean,Interactions2mean,Interactions1Prevmean,Interactions2Prevmean],[stdmean1,stdmean2,stdPrevmean1,stdPrevmean2],'o')
            xticks([1,2,3,4])
            xticklabels({[SITnovelty_mins{SITseshnum,1} 'Trial1'], [SITnovelty_mins{SITseshnum,1} 'Trial2'],[SITnovelty_mins{previousSesh(1),1} 'Trial1'], [SITnovelty_mins{previousSesh(1),1} 'Trial2']})
            xlim([0.8 4.2])
        else
            scatter([1;2;3;4],[Interactions1Prevmean;Interactions2Prevmean;Interactions1mean;Interactions2mean])
            errorbar([1,2,3,4],[Interactions1Prevmean,Interactions2Prevmean,Interactions1mean,Interactions2mean],[stdPrevmean1,stdPrevmean2,stdmean1,stdmean2],'o')
            xticks([1,2,3,4])
            xticklabels({[SITnovelty_mins{previousSesh(1),1} 'Trial1'], [SITnovelty_mins{previousSesh(1),1} 'Trial2'],[SITnovelty_mins{SITseshnum,1} 'Trial1'], [SITnovelty_mins{SITseshnum,1} 'Trial2']})
            xlim([0.8 4.2])
        end
        pass = [];
        if contains(SITnovelty_mins{SITseshnum,1},'Same')
            if SITnovelty_mins{SITseshnum,2} == 30
                diffZonemean1_30(mousecount,1) = Interactions1Prevmean;
                diffZonemean2_30(mousecount,1) = Interactions2Prevmean;
                sameZonemean1_30(mousecount,1) = Interactions1mean;
                sameZonemean2_30(mousecount,1) = Interactions2mean;
                diffZonemean1_30(mousecount,2) = stdPrevmean1;
                diffZonemean2_30(mousecount,2) = stdPrevmean2;
                sameZonemean1_30(mousecount,2) = stdmean1;
                sameZonemean2_30(mousecount,2) = stdmean2;
                
                MeanFiringAllcells30Trial1Same = cat(1,MeanFiringAllcells30Trial1Same,mean(catActivity1,2));
                MeanFiringAllcells30Trial1SameRaw = cat(1,MeanFiringAllcells30Trial1SameRaw,mean(catActivity1Raw,2));
                MeanFiringAllcells30Trial2Same = cat(1,MeanFiringAllcells30Trial2Same,mean(catActivity2,2));
                MeanFiringAllcells30Trial2SameRaw = cat(1,MeanFiringAllcells30Trial2SameRaw,mean(catActivity2Raw,2));
                MeanFiringAllcells30Trial1Diff = cat(1,MeanFiringAllcells30Trial1Diff,mean(prev4b.catActivity1,2));
                MeanFiringAllcells30Trial1DiffRaw = cat(1,MeanFiringAllcells30Trial1DiffRaw,mean(prev4b.catActivity1Raw,2));
                MeanFiringAllcells30Trial2Diff = cat(1,MeanFiringAllcells30Trial2Diff,mean(prev4b.catActivity2,2));
                MeanFiringAllcells30Trial2DiffRaw = cat(1,MeanFiringAllcells30Trial2DiffRaw,mean(prev4b.catActivity2Raw,2));
                
                MeanFiringAllcells30Hab1Same = cat(1,MeanFiringAllcells30Hab1Same,mean(catActivityHab1,2));
                MeanFiringAllcells30Hab1SameRaw = cat(1,MeanFiringAllcells30Hab1SameRaw,mean(catActivityHab1Raw,2));
                MeanFiringAllcells30Hab2Same = cat(1,MeanFiringAllcells30Hab2Same,mean(catActivityHab2,2));
                MeanFiringAllcells30Hab2SameRaw = cat(1,MeanFiringAllcells30Hab2SameRaw,mean(catActivityHab2Raw,2));
                MeanFiringAllcells30Hab1Diff = cat(1,MeanFiringAllcells30Hab1Diff,mean(prev4b.catActivityHab1,2));
                MeanFiringAllcells30Hab1DiffRaw = cat(1,MeanFiringAllcells30Hab1DiffRaw,mean(prev4b.catActivityHab1Raw,2));
                MeanFiringAllcells30Hab2Diff = cat(1,MeanFiringAllcells30Hab2Diff,mean(prev4b.catActivityHab2,2));
                MeanFiringAllcells30Hab2DiffRaw = cat(1,MeanFiringAllcells30Hab2DiffRaw,mean(prev4b.catActivityHab2Raw,2));
                
                MeanFiringAllcells30Trial1Same_passed = cat(1,MeanFiringAllcells30Trial1Same,mean(catActivity1(p1z,:),2));                
                MeanFiringAllcells30Trial2Same_passed = cat(1,MeanFiringAllcells30Trial2Same,mean(catActivity2(p2z,:),2));
                MeanFiringAllcells30Trial1Diff_passed = cat(1,MeanFiringAllcells30Trial1Diff,mean(prev4b.catActivity1(p1prevz,:),2));                
                MeanFiringAllcells30Trial2Diff_passed = cat(1,MeanFiringAllcells30Trial2Diff,mean(prev4b.catActivity2(p2prevz,:),2));                                
                MeanFiringAllcells30Hab1Same_passed = cat(1,MeanFiringAllcells30Hab1Same_passed,mean(catActivityHab1(ph1z,:),2));                
                MeanFiringAllcells30Hab2Same_passed = cat(1,MeanFiringAllcells30Hab2Same_passed,mean(catActivityHab2(ph2z,:),2));                
                MeanFiringAllcells30Hab1Diff_passed = cat(1,MeanFiringAllcells30Hab1Diff_passed,mean(prev4b.catActivityHab1(ph1prevz,:),2));                
                MeanFiringAllcells30Hab2Diff_passed = cat(1,MeanFiringAllcells30Hab2Diff_passed,mean(prev4b.catActivityHab2(ph2prevz,:),2));     
                                                
                %colour coated non-cellreg
                pass = passNonCellReg_same(p1z,ph1z,p2z,ph2z,p1z_not,ph1z_not,p2z_not,ph2z_not);
                %Trial1 Same vs Trial2 Same
                MeanFiring_30_S1S2 = cat(1,MeanFiring_30_S1S2,[mean(catActivity1(pass.pass_S1_S2,:),2),mean(catActivity2(pass.pass_S1_S2,:),2)]);
                MeanFiring_30_Sn1S2 = cat(1,MeanFiring_30_Sn1S2,[mean(catActivity1(pass.pass_Sn1_S2,:),2),mean(catActivity2(pass.pass_Sn1_S2,:),2)]);
                MeanFiring_30_S1Sn2 = cat(1,MeanFiring_30_S1Sn2,[mean(catActivity1(pass.pass_S1_Sn2,:),2),mean(catActivity2(pass.pass_S1_Sn2,:),2)]);
                MeanFiring_30_Sn1Sn2 = cat(1,MeanFiring_30_Sn1Sn2,[mean(catActivity1(pass.pass_Sn1_Sn2,:),2),mean(catActivity2(pass.pass_Sn1_Sn2,:),2)]);
                %Hab1 Same vs Trial2 Same
                MeanFiring_30_Sh1S2 = cat(1,MeanFiring_30_Sh1S2,[mean(catActivityHab1(pass.pass_Sh1_S2,:),2),mean(catActivity2(pass.pass_Sh1_S2,:),2)]);
                MeanFiring_30_Snh1S2 = cat(1,MeanFiring_30_Snh1S2,[mean(catActivityHab1(pass.pass_Snh1_S2,:),2),mean(catActivity2(pass.pass_Snh1_S2,:),2)]);
                MeanFiring_30_Sh1Sn2 = cat(1,MeanFiring_30_Sh1Sn2,[mean(catActivityHab1(pass.pass_Sh1_Sn2,:),2),mean(catActivity2(pass.pass_Sh1_Sn2,:),2)]);
                MeanFiring_30_Snh1Sn2 = cat(1,MeanFiring_30_Snh1Sn2,[mean(catActivityHab1(pass.pass_Snh1_Sn2,:),2),mean(catActivity2(pass.pass_Snh1_Sn2,:),2)]);
                %Trial1 Same vs Hab2 Same
                MeanFiring_30_S1Sh2 = cat(1,MeanFiring_30_S1Sh2,[mean(catActivity1(pass.pass_S1_Sh2,:),2),mean(catActivityHab2(pass.pass_S1_Sh2,:),2)]);
                MeanFiring_30_Sn1Sh2 = cat(1,MeanFiring_30_Sn1Sh2,[mean(catActivity1(pass.pass_Sn1_Sh2,:),2),mean(catActivityHab2(pass.pass_Sn1_Sh2,:),2)]);
                MeanFiring_30_S1Snh2 = cat(1,MeanFiring_30_S1Snh2,[mean(catActivity1(pass.pass_S1_Snh2,:),2),mean(catActivityHab2(pass.pass_S1_Snh2,:),2)]);
                MeanFiring_30_Sn1Snh2 = cat(1,MeanFiring_30_Sn1Snh2,[mean(catActivity1(pass.pass_Sn1_Snh2,:),2),mean(catActivityHab2(pass.pass_Sn1_Snh2,:),2)]);
                %Hab1 Same vs Hab2 Same
                MeanFiring_30_Sh1Sh2 = cat(1,MeanFiring_30_Sh1Sh2,[mean(catActivityHab1(pass.pass_Sh1_Sh2,:),2),mean(catActivityHab2(pass.pass_Sh1_Sh2,:),2)]);
                MeanFiring_30_Snh1Sh2 = cat(1,MeanFiring_30_Snh1Sh2,[mean(catActivityHab1(pass.pass_Snh1_Sh2,:),2),mean(catActivityHab2(pass.pass_Snh1_Sh2,:),2)]);
                MeanFiring_30_Sh1Snh2 = cat(1,MeanFiring_30_Sh1Snh2,[mean(catActivityHab1(pass.pass_Sh1_Snh2,:),2),mean(catActivityHab2(pass.pass_Sh1_Snh2,:),2)]);
                MeanFiring_30_Snh1Snh2 = cat(1,MeanFiring_30_Snh1Snh2,[mean(catActivityHab1(pass.pass_Snh1_Snh2,:),2),mean(catActivityHab2(pass.pass_Snh1_Snh2,:),2)]);
                
                pass = passNonCellReg_diff(p1prevz,ph1prevz,p2prevz,ph2prevz,p1prevz_not,ph1prevz_not,p2prevz_not,ph2prevz_not);
                %Trial1 Diff vs Trial2 Diff
                MeanFiring_30_D1D2 = cat(1,MeanFiringCellReg_30_D1D2,[mean(prev4b.catActivity1(pass.pass_D1_D2,:),2),mean(prev4b.catActivity2(pass.pass_D1_D2,:),2)]);
                MeanFiring_30_Dn1D2 = cat(1,MeanFiringCellReg_30_Dn1D2,[mean(prev4b.catActivity1(pass.pass_Dn1_D2,:),2),mean(prev4b.catActivity2(pass.pass_Dn1_D2,:),2)]);
                MeanFiring_30_D1Dn2 = cat(1,MeanFiringCellReg_30_D1Dn2,[mean(prev4b.catActivity1(pass.pass_D1_Dn2,:),2),mean(prev4b.catActivity2(pass.pass_D1_Dn2,:),2)]);
                MeanFiring_30_Dn1Dn2 = cat(1,MeanFiringCellReg_30_Dn1Dn2,[mean(prev4b.catActivity1(pass.pass_Dn1_Dn2,:),2),mean(prev4b.catActivity2(pass.pass_Dn1_Dn2,:),2)]);
                %Hab1 Diff vs Trial2 Diff
                MeanFiring_30_Dh1D2 = cat(1,MeanFiring_30_Dh1D2,[mean(prev4b.catActivityHab1(pass.pass_Dh1_D2,:),2),mean(prev4b.catActivity2(pass.pass_Dh1_D2,:),2)]);
                MeanFiring_30_Dnh1D2 = cat(1,MeanFiring_30_Dnh1D2,[mean(prev4b.catActivityHab1(pass.pass_Dnh1_D2,:),2),mean(prev4b.catActivity2(pass.pass_Dnh1_D2,:),2)]);
                MeanFiring_30_Dh1Dn2 = cat(1,MeanFiring_30_Dh1Dn2,[mean(prev4b.catActivityHab1(pass.pass_Dh1_Dn2,:),2),mean(prev4b.catActivity2(pass.pass_Dh1_Dn2,:),2)]);
                MeanFiring_30_Dnh1Dn2 = cat(1,MeanFiring_30_Dnh1Dn2,[mean(prev4b.catActivityHab1(pass.pass_Dnh1_Dn2,:),2),mean(prev4b.catActivity2(pass.pass_Dnh1_Dn2,:),2)]);
                %Trial1 Diff vs Hab2 Diff
                MeanFiring_30_D1Dh2 = cat(1,MeanFiring_30_D1Dh2,[mean(prev4b.catActivity1(pass.pass_D1_Dh2,:),2),mean(prev4b.catActivityHab2(pass.pass_D1_Dh2,:),2)]);
                MeanFiring_30_Dn1Dh2 = cat(1,MeanFiring_30_Dn1Dh2,[mean(prev4b.catActivity1(pass.pass_Dn1_Dh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dn1_Dh2,:),2)]);
                MeanFiring_30_D1Dnh2 = cat(1,MeanFiring_30_D1Dnh2,[mean(prev4b.catActivity1(pass.pass_D1_Dnh2,:),2),mean(prev4b.catActivityHab2(pass.pass_D1_Dnh2,:),2)]);
                MeanFiring_30_Dn1Dnh2 = cat(1,MeanFiring_30_Dn1Dnh2,[mean(prev4b.catActivity1(pass.pass_Dn1_Dnh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dn1_Dnh2,:),2)]);
                %Hab1 Diff vs Hab2 Diff
                MeanFiring_30_Dh1Dh2 = cat(1,MeanFiring_30_Dh1Dh2,[mean(prev4b.catActivityHab1(pass.pass_Dh1_Dh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dh1_Dh2,:),2)]);
                MeanFiring_30_Dnh1Dh2 = cat(1,MeanFiring_30_Dnh1Dh2,[mean(prev4b.catActivityHab1(pass.pass_Dnh1_Dh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dnh1_Dh2,:),2)]);
                MeanFiring_30_Dh1Dnh2 = cat(1,MeanFiring_30_Dh1Dnh2,[mean(prev4b.catActivityHab1(pass.pass_Dh1_Dnh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dh1_Dnh2,:),2)]);
                MeanFiring_30_Dnh1Dnh2 = cat(1,MeanFiring_30_Dnh1Dnh2,[mean(prev4b.catActivityHab1(pass.pass_Dnh1_Dnh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dnh1_Dnh2,:),2)]);
                
            else
                diffZonemean1_120(mousecount,1) = Interactions1Prevmean;
                diffZonemean2_120(mousecount,1) = Interactions2Prevmean;
                sameZonemean1_120(mousecount,1) = Interactions1mean;
                sameZonemean2_120(mousecount,1) = Interactions2mean;
                diffZonemean1_120(mousecount,2) = stdPrevmean1;
                diffZonemean2_120(mousecount,2) = stdPrevmean2;
                sameZonemean1_120(mousecount,2) = stdmean1;
                sameZonemean2_120(mousecount,2) = stdmean2;
                
                MeanFiringAllcells120Trial1Same = cat(1,MeanFiringAllcells120Trial1Same,mean(catActivity1,2));
                MeanFiringAllcells120Trial1SameRaw = cat(1,MeanFiringAllcells120Trial1SameRaw,mean(catActivity1Raw,2));
                MeanFiringAllcells120Trial2Same = cat(1,MeanFiringAllcells120Trial2Same,mean(catActivity2,2));
                MeanFiringAllcells120Trial2SameRaw = cat(1,MeanFiringAllcells120Trial2SameRaw,mean(catActivity2Raw,2));
                MeanFiringAllcells120Trial1Diff = cat(1,MeanFiringAllcells120Trial1Diff,mean(prev4b.catActivity1,2));
                MeanFiringAllcells120Trial1DiffRaw = cat(1,MeanFiringAllcells120Trial1DiffRaw,mean(prev4b.catActivity1Raw,2));
                MeanFiringAllcells120Trial2Diff = cat(1,MeanFiringAllcells120Trial2Diff,mean(prev4b.catActivity2,2));
                MeanFiringAllcells120Trial2DiffRaw = cat(1,MeanFiringAllcells120Trial2DiffRaw,mean(prev4b.catActivity2Raw,2));
                
                MeanFiringAllcells120Hab1Same = cat(1,MeanFiringAllcells120Hab1Same,mean(catActivityHab1,2));
                MeanFiringAllcells120Hab1SameRaw = cat(1,MeanFiringAllcells120Hab1SameRaw,mean(catActivityHab1Raw,2));
                MeanFiringAllcells120Hab2Same = cat(1,MeanFiringAllcells120Hab2Same,mean(catActivityHab2,2));
                MeanFiringAllcells120Hab2SameRaw = cat(1,MeanFiringAllcells120Hab2SameRaw,mean(catActivityHab2Raw,2));
                MeanFiringAllcells120Hab1Diff = cat(1,MeanFiringAllcells120Hab1Diff,mean(prev4b.catActivityHab1,2));
                MeanFiringAllcells120Hab1DiffRaw = cat(1,MeanFiringAllcells120Hab1DiffRaw,mean(prev4b.catActivityHab1Raw,2));
                MeanFiringAllcells120Hab2Diff = cat(1,MeanFiringAllcells120Hab2Diff,mean(prev4b.catActivityHab2,2));
                MeanFiringAllcells120Hab2DiffRaw = cat(1,MeanFiringAllcells120Hab2DiffRaw,mean(prev4b.catActivityHab2Raw,2));
                
                MeanFiringAllcells120Trial1Same_passed = cat(1,MeanFiringAllcells120Trial1Same_passed,mean(catActivity1(p1z,:),2));                
                MeanFiringAllcells120Trial2Same_passed = cat(1,MeanFiringAllcells120Trial2Same_passed,mean(catActivity2(p2z,:),2));                
                MeanFiringAllcells120Trial1Diff_passed = cat(1,MeanFiringAllcells120Trial1Diff_passed,mean(prev4b.catActivity1(p1prevz,:),2));                
                MeanFiringAllcells120Trial2Diff_passed = cat(1,MeanFiringAllcells120Trial2Diff_passed,mean(prev4b.catActivity2(p2prevz,:),2));                                
                MeanFiringAllcells120Hab1Same_passed = cat(1,MeanFiringAllcells120Hab1Same_passed,mean(catActivityHab1(ph1z,:),2));                
                MeanFiringAllcells120Hab2Same_passed = cat(1,MeanFiringAllcells120Hab2Same_passed,mean(catActivityHab2(ph2z,:),2));                
                MeanFiringAllcells120Hab1Diff_passed = cat(1,MeanFiringAllcells120Hab1Diff_passed,mean(prev4b.catActivityHab1(ph1prevz,:),2));                
                MeanFiringAllcells120Hab2Diff_passed = cat(1,MeanFiringAllcells120Hab2Diff_passed,mean(prev4b.catActivityHab2(ph2prevz,:),2));
                
               %colour coated non-cellreg
                pass = passNonCellReg_same(p1z,ph1z,p2z,ph2z,p1z_not,ph1z_not,p2z_not,ph2z_not);
                %Trial1 Same vs Trial2 Same
                MeanFiring_120_S1S2 = cat(1,MeanFiring_120_S1S2,[mean(catActivity1(pass.pass_S1_S2,:),2),mean(catActivity2(pass.pass_S1_S2,:),2)]);
                MeanFiring_120_Sn1S2 = cat(1,MeanFiring_120_Sn1S2,[mean(catActivity1(pass.pass_Sn1_S2,:),2),mean(catActivity2(pass.pass_Sn1_S2,:),2)]);
                MeanFiring_120_S1Sn2 = cat(1,MeanFiring_120_S1Sn2,[mean(catActivity1(pass.pass_S1_Sn2,:),2),mean(catActivity2(pass.pass_S1_Sn2,:),2)]);
                MeanFiring_120_Sn1Sn2 = cat(1,MeanFiring_120_Sn1Sn2,[mean(catActivity1(pass.pass_Sn1_Sn2,:),2),mean(catActivity2(pass.pass_Sn1_Sn2,:),2)]);
                %Hab1 Same vs Trial2 Same
                MeanFiring_120_Sh1S2 = cat(1,MeanFiring_120_Sh1S2,[mean(catActivityHab1(pass.pass_Sh1_S2,:),2),mean(catActivity2(pass.pass_Sh1_S2,:),2)]);
                MeanFiring_120_Snh1S2 = cat(1,MeanFiring_120_Snh1S2,[mean(catActivityHab1(pass.pass_Snh1_S2,:),2),mean(catActivity2(pass.pass_Snh1_S2,:),2)]);
                MeanFiring_120_Sh1Sn2 = cat(1,MeanFiring_120_Sh1Sn2,[mean(catActivityHab1(pass.pass_Sh1_Sn2,:),2),mean(catActivity2(pass.pass_Sh1_Sn2,:),2)]);
                MeanFiring_120_Snh1Sn2 = cat(1,MeanFiring_120_Snh1Sn2,[mean(catActivityHab1(pass.pass_Snh1_Sn2,:),2),mean(catActivity2(pass.pass_Snh1_Sn2,:),2)]);
                %Trial1 Same vs Hab2 Same
                MeanFiring_120_S1Sh2 = cat(1,MeanFiring_120_S1Sh2,[mean(catActivity1(pass.pass_S1_Sh2(:,1),:),2),mean(catActivityHab2(pass.pass_S1_Sh2,:),2)]);
                MeanFiring_120_Sn1Sh2 = cat(1,MeanFiring_120_Sn1Sh2,[mean(catActivity1(pass.pass_Sn1_Sh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sn1_Sh2,:),2)]);
                MeanFiring_120_S1Snh2 = cat(1,MeanFiring_120_S1Snh2,[mean(catActivity1(pass.pass_S1_Snh2(:,1),:),2),mean(catActivityHab2(pass.pass_S1_Snh2,:),2)]);
                MeanFiring_120_Sn1Snh2 = cat(1,MeanFiring_120_Sn1Snh2,[mean(catActivity1(pass.pass_Sn1_Snh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sn1_Snh2,:),2)]);
                %Hab1 Same vs Hab2 Same
                MeanFiring_120_Sh1Sh2 = cat(1,MeanFiring_120_Sh1Sh2,[mean(catActivityHab1(pass.pass_Sh1_Sh2,:),2),mean(catActivityHab2(pass.pass_Sh1_Sh2,:),2)]);
                MeanFiring_120_Snh1Sh2 = cat(1,MeanFiring_120_Snh1Sh2,[mean(catActivityHab1(pass.pass_Snh1_Sh2,:),2),mean(catActivityHab2(pass.pass_Snh1_Sh2,:),2)]);
                MeanFiring_120_Sh1Snh2 = cat(1,MeanFiring_120_Sh1Snh2,[mean(catActivityHab1(pass.pass_Sh1_Snh2,:),2),mean(catActivityHab2(pass.pass_Sh1_Snh2,:),2)]);
                MeanFiring_120_Snh1Snh2 = cat(1,MeanFiring_120_Snh1Snh2,[mean(catActivityHab1(pass.pass_Snh1_Snh2,:),2),mean(catActivityHab2(pass.pass_Snh1_Snh2,:),2)]);
                
                pass = passNonCellReg_diff(p1prevz,ph1prevz,p2prevz,ph2prevz,p1prevz_not,ph1prevz_not,p2prevz_not,ph2prevz_not);
                %Trial1 Diff vs Trial2 Diff
                MeanFiring_120_D1D2 = cat(1,MeanFiringCellReg_30_D1D2,[mean(prev4b.catActivity1(pass.pass_D1_D2,:),2),mean(prev4b.catActivity2(pass.pass_D1_D2,:),2)]);
                MeanFiring_120_Dn1D2 = cat(1,MeanFiringCellReg_30_Dn1D2,[mean(prev4b.catActivity1(pass.pass_Dn1_D2,:),2),mean(prev4b.catActivity2(pass.pass_Dn1_D2,:),2)]);
                MeanFiring_120_D1Dn2 = cat(1,MeanFiringCellReg_30_D1Dn2,[mean(prev4b.catActivity1(pass.pass_D1_Dn2,:),2),mean(prev4b.catActivity2(pass.pass_D1_Dn2,:),2)]);
                MeanFiring_120_Dn1Dn2 = cat(1,MeanFiringCellReg_30_Dn1Dn2,[mean(prev4b.catActivity1(pass.pass_Dn1_Dn2,:),2),mean(prev4b.catActivity2(pass.pass_Dn1_Dn2,:),2)]);
                %Hab1 Diff vs Trial2 Diff
                MeanFiring_120_Dh1D2 = cat(1,MeanFiring_120_Dh1D2,[mean(prev4b.catActivityHab1(pass.pass_Dh1_D2,:),2),mean(prev4b.catActivity2(pass.pass_Dh1_D2,:),2)]);
                MeanFiring_120_Dnh1D2 = cat(1,MeanFiring_120_Dnh1D2,[mean(prev4b.catActivityHab1(pass.pass_Dnh1_D2,:),2),mean(prev4b.catActivity2(pass.pass_Dnh1_D2,:),2)]);
                MeanFiring_120_Dh1Dn2 = cat(1,MeanFiring_120_Dh1Dn2,[mean(prev4b.catActivityHab1(pass.pass_Dh1_Dn2,:),2),mean(prev4b.catActivity2(pass.pass_Dh1_Dn2,:),2)]);
                MeanFiring_120_Dnh1Dn2 = cat(1,MeanFiring_120_Dnh1Dn2,[mean(prev4b.catActivityHab1(pass.pass_Dnh1_Dn2,:),2),mean(prev4b.catActivity2(pass.pass_Dnh1_Dn2,:),2)]);
                %Trial1 Diff vs Hab2 Diff
                MeanFiring_120_D1Dh2 = cat(1,MeanFiring_120_D1Dh2,[mean(prev4b.catActivity1(pass.pass_D1_Dh2,:),2),mean(prev4b.catActivityHab2(pass.pass_D1_Dh2,:),2)]);
                MeanFiring_120_Dn1Dh2 = cat(1,MeanFiring_120_Dn1Dh2,[mean(prev4b.catActivity1(pass.pass_Dn1_Dh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dn1_Dh2,:),2)]);
                MeanFiring_120_D1Dnh2 = cat(1,MeanFiring_120_D1Dnh2,[mean(prev4b.catActivity1(pass.pass_D1_Dnh2,:),2),mean(prev4b.catActivityHab2(pass.pass_D1_Dnh2,:),2)]);
                MeanFiring_120_Dn1Dnh2 = cat(1,MeanFiring_120_Dn1Dnh2,[mean(prev4b.catActivity1(pass.pass_Dn1_Dnh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dn1_Dnh2,:),2)]);
                %Hab1 Diff vs Hab2 Diff
                MeanFiring_120_Dh1Dh2 = cat(1,MeanFiring_120_Dh1Dh2,[mean(prev4b.catActivityHab1(pass.pass_Dh1_Dh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dh1_Dh2,:),2)]);
                MeanFiring_120_Dnh1Dh2 = cat(1,MeanFiring_120_Dnh1Dh2,[mean(prev4b.catActivityHab1(pass.pass_Dnh1_Dh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dnh1_Dh2,:),2)]);
                MeanFiring_120_Dh1Dnh2 = cat(1,MeanFiring_120_Dh1Dnh2,[mean(prev4b.catActivityHab1(pass.pass_Dh1_Dnh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dh1_Dnh2,:),2)]);
                MeanFiring_120_Dnh1Dnh2 = cat(1,MeanFiring_120_Dnh1Dnh2,[mean(prev4b.catActivityHab1(pass.pass_Dnh1_Dnh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Dnh1_Dnh2,:),2)]);
                
            end
            
        else
            if SITnovelty_mins{SITseshnum,2} == 30
                sameZonemean1_30(mousecount,1) = Interactions1Prevmean;
                sameZonemean2_30(mousecount,1) = Interactions2Prevmean;
                diffZonemean1_30(mousecount,1) = Interactions1mean;
                diffZonemean2_30(mousecount,1) = Interactions2mean;
                sameZonemean1_30(mousecount,2) = stdPrevmean1;
                sameZonemean2_30(mousecount,2) = stdPrevmean2;
                diffZonemean1_30(mousecount,2) = stdmean1;
                diffZonemean2_30(mousecount,2) = stdmean2;
                
                MeanFiringAllcells30Trial1Same = cat(1,MeanFiringAllcells30Trial1Same,mean(prev4b.catActivity1,2));
                MeanFiringAllcells30Trial1SameRaw = cat(1,MeanFiringAllcells30Trial1SameRaw,mean(prev4b.catActivity1Raw,2));
                MeanFiringAllcells30Trial2Same = cat(1,MeanFiringAllcells30Trial2Same,mean(prev4b.catActivity2,2));
                MeanFiringAllcells30Trial2SameRaw = cat(1,MeanFiringAllcells30Trial2SameRaw,mean(prev4b.catActivity2Raw,2));
                MeanFiringAllcells30Trial1Diff = cat(1,MeanFiringAllcells30Trial1Diff,mean(catActivity1,2));
                MeanFiringAllcells30Trial1DiffRaw = cat(1,MeanFiringAllcells30Trial1DiffRaw,mean(catActivity1Raw,2));
                MeanFiringAllcells30Trial2Diff = cat(1,MeanFiringAllcells30Trial2Diff,mean(catActivity2,2));
                MeanFiringAllcells30Trial2DiffRaw = cat(1,MeanFiringAllcells30Trial2DiffRaw,mean(catActivity2Raw,2));
                
                MeanFiringAllcells30Hab1Same = cat(1,MeanFiringAllcells30Hab1Same,mean(prev4b.catActivityHab1,2));
                MeanFiringAllcells30Hab1SameRaw = cat(1,MeanFiringAllcells30Hab1SameRaw,mean(prev4b.catActivityHab1Raw,2));
                MeanFiringAllcells30Hab2Same = cat(1,MeanFiringAllcells30Hab2Same,mean(prev4b.catActivityHab2,2));
                MeanFiringAllcells30Hab2SameRaw = cat(1,MeanFiringAllcells30Hab2SameRaw,mean(prev4b.catActivityHab2Raw,2));
                MeanFiringAllcells30Hab1Diff = cat(1,MeanFiringAllcells30Hab1Diff,mean(catActivityHab1,2));
                MeanFiringAllcells30Hab1DiffRaw = cat(1,MeanFiringAllcells30Hab1DiffRaw,mean(catActivityHab1Raw,2));
                MeanFiringAllcells30Hab2Diff = cat(1,MeanFiringAllcells30Hab2Diff,mean(catActivityHab2,2));
                MeanFiringAllcells30Hab2DiffRaw = cat(1,MeanFiringAllcells30Hab2DiffRaw,mean(catActivityHab2Raw,2));
                
                MeanFiringAllcells30Trial1Same_passed = cat(1,MeanFiringAllcells30Trial1Same_passed,mean(prev4b.catActivity1(p1prevz,:),2));               
                MeanFiringAllcells30Trial2Same_passed = cat(1,MeanFiringAllcells30Trial2Same_passed,mean(prev4b.catActivity2(p2prevz,:),2));                
                MeanFiringAllcells30Trial1Diff_passed = cat(1,MeanFiringAllcells30Trial1Diff_passed,mean(catActivity1(p1z,:),2));               
                MeanFiringAllcells30Trial2Diff_passed = cat(1,MeanFiringAllcells30Trial2Diff_passed,mean(catActivity2(p2z,:),2));                                
                MeanFiringAllcells30Hab1Same_passed = cat(1,MeanFiringAllcells30Hab1Same_passed,mean(prev4b.catActivityHab1(ph1prevz,:),2));                
                MeanFiringAllcells30Hab2Same_passed = cat(1,MeanFiringAllcells30Hab2Same_passed,mean(prev4b.catActivityHab2(ph2prevz,:),2));                
                MeanFiringAllcells30Hab1Diff_passed = cat(1,MeanFiringAllcells30Hab1Diff_passed,mean(catActivityHab1(ph1z,:),2));                
                MeanFiringAllcells30Hab2Diff_passed = cat(1,MeanFiringAllcells30Hab2Diff_passed,mean(catActivityHab2(ph2z,:),2));  
                
                %colour coated non-cellreg
                pass = passNonCellReg_same(p1prevz,ph1prevz,p2prevz,ph2prevz,p1prevz_not,ph1prevz_not,p2prevz_not,ph2prevz_not);
                %Trial1 Same vs Trial2 Same
                MeanFiring_30_S1S2 = cat(1,MeanFiring_30_S1S2,[mean(prev4b.catActivity1(pass.pass_S1_S2,:),2),mean(prev4b.catActivity2(pass.pass_S1_S2,:),2)]);
                MeanFiring_30_Sn1S2 = cat(1,MeanFiring_30_Sn1S2,[mean(prev4b.catActivity1(pass.pass_Sn1_S2,:),2),mean(prev4b.catActivity2(pass.pass_Sn1_S2,:),2)]);
                MeanFiring_30_S1Sn2 = cat(1,MeanFiring_30_S1Sn2,[mean(prev4b.catActivity1(pass.pass_S1_Sn2,:),2),mean(prev4b.catActivity2(pass.pass_S1_Sn2,:),2)]);
                MeanFiring_30_Sn1Sn2 = cat(1,MeanFiring_30_Sn1Sn2,[mean(prev4b.catActivity1(pass.pass_Sn1_Sn2,:),2),mean(prev4b.catActivity2(pass.pass_Sn1_Sn2,:),2)]);
                %Hab1 Same vs Trial2 Same
                MeanFiring_30_Sh1S2 = cat(1,MeanFiring_30_Sh1S2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_S2,:),2),mean(prev4b.catActivity2(pass.pass_Sh1_S2,:),2)]);
                MeanFiring_30_Snh1S2 = cat(1,MeanFiring_30_Snh1S2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_S2,:),2),mean(prev4b.catActivity2(pass.pass_Snh1_S2,:),2)]);
                MeanFiring_30_Sh1Sn2 = cat(1,MeanFiring_30_Sh1Sn2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Sn2,:),2),mean(prev4b.catActivity2(pass.pass_Sh1_Sn2,:),2)]);
                MeanFiring_30_Snh1Sn2 = cat(1,MeanFiring_30_Snh1Sn2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Sn2,:),2),mean(prev4b.catActivity2(pass.pass_Snh1_Sn2,:),2)]);
                %Trial1 Same vs Hab2 Same
                MeanFiring_30_S1Sh2 = cat(1,MeanFiring_30_S1Sh2,[mean(prev4b.catActivity1(pass.pass_S1_Sh2,:),2),mean(prev4b.catActivityHab2(pass.pass_S1_Sh2,:),2)]);
                MeanFiring_30_Sn1Sh2 = cat(1,MeanFiring_30_Sn1Sh2,[mean(prev4b.catActivity1(pass.pass_Sn1_Sh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Sn1_Sh2,:),2)]);
                MeanFiring_30_S1Snh2 = cat(1,MeanFiring_30_S1Snh2,[mean(prev4b.catActivity1(pass.pass_S1_Snh2,:),2),mean(prev4b.catActivityHab2(pass.pass_S1_Snh2,:),2)]);
                MeanFiring_30_Sn1Snh2 = cat(1,MeanFiring_30_Sn1Snh2,[mean(prev4b.catActivity1(pass.pass_Sn1_Snh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Sn1_Snh2,:),2)]);
                %Hab1 Same vs Hab2 Same
                MeanFiring_30_Sh1Sh2 = cat(1,MeanFiring_30_Sh1Sh2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Sh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Sh1_Sh2,:),2)]);
                MeanFiring_30_Snh1Sh2 = cat(1,MeanFiring_30_Snh1Sh2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Sh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Snh1_Sh2,:),2)]);
                MeanFiring_30_Sh1Snh2 = cat(1,MeanFiring_30_Sh1Snh2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Snh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Sh1_Snh2,:),2)]);
                MeanFiring_30_Snh1Snh2 = cat(1,MeanFiring_30_Snh1Snh2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Snh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Snh1_Snh2,:),2)]);
                
                pass = passNonCellReg_diff(p1z,ph1z,p2z,ph2z,p1z_not,ph1z_not,p2z_not,ph2z_not);
                %Trial1 Diff vs Trial2 Diff
                MeanFiring_30_D1D2 = cat(1,MeanFiringCellReg_30_D1D2,[mean(catActivity1(pass.pass_D1_D2,:),2),mean(catActivity2(pass.pass_D1_D2,:),2)]);
                MeanFiring_30_Dn1D2 = cat(1,MeanFiringCellReg_30_Dn1D2,[mean(catActivity1(pass.pass_Dn1_D2,:),2),mean(catActivity2(pass.pass_Dn1_D2,:),2)]);
                MeanFiring_30_D1Dn2 = cat(1,MeanFiringCellReg_30_D1Dn2,[mean(catActivity1(pass.pass_D1_Dn2,:),2),mean(catActivity2(pass.pass_D1_Dn2,:),2)]);
                MeanFiring_30_Dn1Dn2 = cat(1,MeanFiringCellReg_30_Dn1Dn2,[mean(catActivity1(pass.pass_Dn1_Dn2,:),2),mean(catActivity2(pass.pass_Dn1_Dn2,:),2)]);
                %Hab1 Diff vs Trial2 Diff
                MeanFiring_30_Dh1D2 = cat(1,MeanFiring_30_Dh1D2,[mean(catActivityHab1(pass.pass_Dh1_D2,:),2),mean(catActivity2(pass.pass_Dh1_D2,:),2)]);
                MeanFiring_30_Dnh1D2 = cat(1,MeanFiring_30_Dnh1D2,[mean(catActivityHab1(pass.pass_Dnh1_D2,:),2),mean(catActivity2(pass.pass_Dnh1_D2,:),2)]);
                MeanFiring_30_Dh1Dn2 = cat(1,MeanFiring_30_Dh1Dn2,[mean(catActivityHab1(pass.pass_Dh1_Dn2,:),2),mean(catActivity2(pass.pass_Dh1_Dn2,:),2)]);
                MeanFiring_30_Dnh1Dn2 = cat(1,MeanFiring_30_Dnh1Dn2,[mean(catActivityHab1(pass.pass_Dnh1_Dn2,:),2),mean(catActivity2(pass.pass_Dnh1_Dn2,:),2)]);
                %Trial1 Diff vs Hab2 Diff
                MeanFiring_30_D1Dh2 = cat(1,MeanFiring_30_D1Dh2,[mean(catActivity1(pass.pass_D1_Dh2,:),2),mean(catActivityHab2(pass.pass_D1_Dh2,:),2)]);
                MeanFiring_30_Dn1Dh2 = cat(1,MeanFiring_30_Dn1Dh2,[mean(catActivity1(pass.pass_Dn1_Dh2,:),2),mean(catActivityHab2(pass.pass_Dn1_Dh2,:),2)]);
                MeanFiring_30_D1Dnh2 = cat(1,MeanFiring_30_D1Dnh2,[mean(catActivity1(pass.pass_D1_Dnh2,:),2),mean(catActivityHab2(pass.pass_D1_Dnh2,:),2)]);
                MeanFiring_30_Dn1Dnh2 = cat(1,MeanFiring_30_Dn1Dnh2,[mean(catActivity1(pass.pass_Dn1_Dnh2,:),2),mean(catActivityHab2(pass.pass_Dn1_Dnh2,:),2)]);
                %Hab1 Diff vs Hab2 Diff
                MeanFiring_30_Dh1Dh2 = cat(1,MeanFiring_30_Dh1Dh2,[mean(catActivityHab1(pass.pass_Dh1_Dh2,:),2),mean(catActivityHab2(pass.pass_Dh1_Dh2,:),2)]);
                MeanFiring_30_Dnh1Dh2 = cat(1,MeanFiring_30_Dnh1Dh2,[mean(catActivityHab1(pass.pass_Dnh1_Dh2,:),2),mean(catActivityHab2(pass.pass_Dnh1_Dh2,:),2)]);
                MeanFiring_30_Dh1Dnh2 = cat(1,MeanFiring_30_Dh1Dnh2,[mean(catActivityHab1(pass.pass_Dh1_Dnh2,:),2),mean(catActivityHab2(pass.pass_Dh1_Dnh2,:),2)]);
                MeanFiring_30_Dnh1Dnh2 = cat(1,MeanFiring_30_Dnh1Dnh2,[mean(catActivityHab1(pass.pass_Dnh1_Dnh2,:),2),mean(catActivityHab2(pass.pass_Dnh1_Dnh2,:),2)]);
                                
            else
                sameZonemean1_120(mousecount,1) = Interactions1Prevmean;
                sameZonemean2_120(mousecount,1) = Interactions2Prevmean;
                diffZonemean1_120(mousecount,1) = Interactions1mean;
                diffZonemean2_120(mousecount,1) = Interactions2mean;
                sameZonemean1_120(mousecount,2) = stdPrevmean1;
                sameZonemean2_120(mousecount,2) = stdPrevmean2;
                diffZonemean1_120(mousecount,2) = stdmean1;
                diffZonemean2_120(mousecount,2) = stdmean2;
                
                MeanFiringAllcells120Trial1Same = cat(1,MeanFiringAllcells120Trial1Same,mean(prev4b.catActivity1,2));
                MeanFiringAllcells120Trial1SameRaw = cat(1,MeanFiringAllcells120Trial1SameRaw,mean(prev4b.catActivity1Raw,2));
                MeanFiringAllcells120Trial2Same = cat(1,MeanFiringAllcells120Trial2Same,mean(prev4b.catActivity2,2));
                MeanFiringAllcells120Trial2SameRaw = cat(1,MeanFiringAllcells120Trial2SameRaw,mean(prev4b.catActivity2Raw,2));
                MeanFiringAllcells120Trial1Diff = cat(1,MeanFiringAllcells120Trial1Diff,mean(catActivity1,2));
                MeanFiringAllcells120Trial1DiffRaw = cat(1,MeanFiringAllcells120Trial1DiffRaw,mean(catActivity1Raw,2));
                MeanFiringAllcells120Trial2Diff = cat(1,MeanFiringAllcells120Trial2Diff,mean(catActivity2,2));
                MeanFiringAllcells120Trial2DiffRaw = cat(1,MeanFiringAllcells120Trial2DiffRaw,mean(catActivity2Raw,2));
                
                MeanFiringAllcells120Hab1Same = cat(1,MeanFiringAllcells120Hab1Same,mean(prev4b.catActivityHab1,2));
                MeanFiringAllcells120Hab1SameRaw = cat(1,MeanFiringAllcells120Hab1SameRaw,mean(prev4b.catActivityHab1Raw,2));
                MeanFiringAllcells120Hab2Same = cat(1,MeanFiringAllcells120Hab2Same,mean(prev4b.catActivityHab2,2));
                MeanFiringAllcells120Hab2SameRaw = cat(1,MeanFiringAllcells120Hab2SameRaw,mean(prev4b.catActivityHab2Raw,2));
                MeanFiringAllcells120Hab1Diff = cat(1,MeanFiringAllcells120Hab1Diff,mean(catActivityHab1,2));
                MeanFiringAllcells120Hab1DiffRaw = cat(1,MeanFiringAllcells120Hab1DiffRaw,mean(catActivityHab1Raw,2));
                MeanFiringAllcells120Hab2Diff = cat(1,MeanFiringAllcells120Hab2Diff,mean(catActivityHab2,2));
                MeanFiringAllcells120Hab2DiffRaw = cat(1,MeanFiringAllcells120Hab2DiffRaw,mean(catActivityHab2Raw,2));
                
                MeanFiringAllcells120Trial1Same_passed = cat(1,MeanFiringAllcells120Trial1Same_passed,mean(prev4b.catActivity1(p1prevz,:),2));
                MeanFiringAllcells120Trial2Same_passed = cat(1,MeanFiringAllcells120Trial2Same_passed,mean(prev4b.catActivity2(p2prevz,:),2));
                MeanFiringAllcells120Trial1Diff_passed = cat(1,MeanFiringAllcells120Trial1Diff_passed,mean(catActivity1(p1z,:),2));
                MeanFiringAllcells120Trial2Diff_passed = cat(1,MeanFiringAllcells120Trial2Diff_passed,mean(catActivity2(p2z,:),2));
                MeanFiringAllcells120Hab1Same_passed = cat(1,MeanFiringAllcells120Hab1Same_passed,mean(prev4b.catActivityHab1(ph1prevz,:),2));
                MeanFiringAllcells120Hab2Same_passed = cat(1,MeanFiringAllcells120Hab2Same_passed,mean(prev4b.catActivityHab2(ph2prevz,:),2));
                MeanFiringAllcells120Hab1Diff_passed = cat(1,MeanFiringAllcells120Hab1Diff_passed,mean(catActivityHab1(ph1z,:),2));
                MeanFiringAllcells120Hab2Diff_passed = cat(1,MeanFiringAllcells120Hab2Diff_passed,mean(catActivityHab2(ph2z,:),2));     
                
                %colour coated non-cellreg
                pass = passNonCellReg_same(p1prevz,ph1prevz,p2prevz,ph2prevz,p1prevz_not,ph1prevz_not,p2prevz_not,ph2prevz_not);
                %Trial1 Same vs Trial2 Same
                MeanFiring_120_S1S2 = cat(1,MeanFiring_120_S1S2,[mean(prev4b.catActivity1(pass.pass_S1_S2,:),2),mean(prev4b.catActivity2(pass.pass_S1_S2,:),2)]);
                MeanFiring_120_Sn1S2 = cat(1,MeanFiring_120_Sn1S2,[mean(prev4b.catActivity1(pass.pass_Sn1_S2,:),2),mean(prev4b.catActivity2(pass.pass_Sn1_S2,:),2)]);
                MeanFiring_120_S1Sn2 = cat(1,MeanFiring_120_S1Sn2,[mean(prev4b.catActivity1(pass.pass_S1_Sn2,:),2),mean(prev4b.catActivity2(pass.pass_S1_Sn2,:),2)]);
                MeanFiring_120_Sn1Sn2 = cat(1,MeanFiring_120_Sn1Sn2,[mean(prev4b.catActivity1(pass.pass_Sn1_Sn2,:),2),mean(prev4b.catActivity2(pass.pass_Sn1_Sn2,:),2)]);
                %Hab1 Same vs Trial2 Same
                MeanFiring_120_Sh1S2 = cat(1,MeanFiring_120_Sh1S2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_S2,:),2),mean(prev4b.catActivity2(pass.pass_Sh1_S2,:),2)]);
                MeanFiring_120_Snh1S2 = cat(1,MeanFiring_120_Snh1S2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_S2,:),2),mean(prev4b.catActivity2(pass.pass_Snh1_S2,:),2)]);
                MeanFiring_120_Sh1Sn2 = cat(1,MeanFiring_120_Sh1Sn2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Sn2,:),2),mean(prev4b.catActivity2(pass.pass_Sh1_Sn2,:),2)]);
                MeanFiring_120_Snh1Sn2 = cat(1,MeanFiring_120_Snh1Sn2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Sn2,:),2),mean(prev4b.catActivity2(pass.pass_Snh1_Sn2,:),2)]);
                %Trial1 Same vs Hab2 Same
                MeanFiring_120_S1Sh2 = cat(1,MeanFiring_120_S1Sh2,[mean(prev4b.catActivity1(pass.pass_S1_Sh2,:),2),mean(prev4b.catActivityHab2(pass.pass_S1_Sh2,:),2)]);
                MeanFiring_120_Sn1Sh2 = cat(1,MeanFiring_120_Sn1Sh2,[mean(prev4b.catActivity1(pass.pass_Sn1_Sh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Sn1_Sh2,:),2)]);
                MeanFiring_120_S1Snh2 = cat(1,MeanFiring_120_S1Snh2,[mean(prev4b.catActivity1(pass.pass_S1_Snh2,:),2),mean(prev4b.catActivityHab2(pass.pass_S1_Snh2,:),2)]);
                MeanFiring_120_Sn1Snh2 = cat(1,MeanFiring_120_Sn1Snh2,[mean(prev4b.catActivity1(pass.pass_Sn1_Snh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Sn1_Snh2,:),2)]);
                %Hab1 Same vs Hab2 Same
                MeanFiring_120_Sh1Sh2 = cat(1,MeanFiring_120_Sh1Sh2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Sh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Sh1_Sh2,:),2)]);
                MeanFiring_120_Snh1Sh2 = cat(1,MeanFiring_120_Snh1Sh2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Sh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Snh1_Sh2,:),2)]);
                MeanFiring_120_Sh1Snh2 = cat(1,MeanFiring_120_Sh1Snh2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Snh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Sh1_Snh2,:),2)]);
                MeanFiring_120_Snh1Snh2 = cat(1,MeanFiring_120_Snh1Snh2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Snh2,:),2),mean(prev4b.catActivityHab2(pass.pass_Snh1_Snh2,:),2)]);
                
                pass = passNonCellReg_diff(p1z,ph1z,p2z,ph2z,p1z_not,ph1z_not,p2z_not,ph2z_not);
                %Trial1 Diff vs Trial2 Diff
                MeanFiring_120_D1D2 = cat(1,MeanFiringCellReg_120_D1D2,[mean(catActivity1(pass.pass_D1_D2,:),2),mean(catActivity2(pass.pass_D1_D2,:),2)]);
                MeanFiring_120_Dn1D2 = cat(1,MeanFiringCellReg_120_Dn1D2,[mean(catActivity1(pass.pass_Dn1_D2,:),2),mean(catActivity2(pass.pass_Dn1_D2,:),2)]);
                MeanFiring_120_D1Dn2 = cat(1,MeanFiringCellReg_120_D1Dn2,[mean(catActivity1(pass.pass_D1_Dn2,:),2),mean(catActivity2(pass.pass_D1_Dn2,:),2)]);
                MeanFiring_120_Dn1Dn2 = cat(1,MeanFiringCellReg_120_Dn1Dn2,[mean(catActivity1(pass.pass_Dn1_Dn2,:),2),mean(catActivity2(pass.pass_Dn1_Dn2,:),2)]);
                %Hab1 Diff 1 vs Trial2 Diff
                MeanFiring_120_Dh1D2 = cat(1,MeanFiring_120_Dh1D2,[mean(catActivityHab1(pass.pass_Dh1_D2,:),2),mean(catActivity2(pass.pass_Dh1_D2,:),2)]);
                MeanFiring_120_Dnh1D2 = cat(1,MeanFiring_120_Dnh1D2,[mean(catActivityHab1(pass.pass_Dnh1_D2,:),2),mean(catActivity2(pass.pass_Dnh1_D2,:),2)]);
                MeanFiring_120_Dh1Dn2 = cat(1,MeanFiring_120_Dh1Dn2,[mean(catActivityHab1(pass.pass_Dh1_Dn2,:),2),mean(catActivity2(pass.pass_Dh1_Dn2,:),2)]);
                MeanFiring_120_Dnh1Dn2 = cat(1,MeanFiring_120_Dnh1Dn2,[mean(catActivityHab1(pass.pass_Dnh1_Dn2,:),2),mean(catActivity2(pass.pass_Dnh1_Dn2,:),2)]);
                %Trial1 Diff vs Hab2 Diff
                MeanFiring_120_D1Dh2 = cat(1,MeanFiring_120_D1Dh2,[mean(catActivity1(pass.pass_D1_Dh2,:),2),mean(catActivityHab2(pass.pass_D1_Dh2,:),2)]);
                MeanFiring_120_Dn1Dh2 = cat(1,MeanFiring_120_Dn1Dh2,[mean(catActivity1(pass.pass_Dn1_Dh2,:),2),mean(catActivityHab2(pass.pass_Dn1_Dh2,:),2)]);
                MeanFiring_120_D1Dnh2 = cat(1,MeanFiring_120_D1Dnh2,[mean(catActivity1(pass.pass_D1_Dnh2,:),2),mean(catActivityHab2(pass.pass_D1_Dnh2,:),2)]);
                MeanFiring_120_Dn1Dnh2 = cat(1,MeanFiring_120_Dn1Dnh2,[mean(catActivity1(pass.pass_Dn1_Dnh2,:),2),mean(catActivityHab2(pass.pass_Dn1_Dnh2,:),2)]);
                %Hab1 Diff vs Hab2 Diff
                MeanFiring_120_Dh1Dh2 = cat(1,MeanFiring_120_Dh1Dh2,[mean(catActivityHab1(pass.pass_Dh1_Dh2,:),2),mean(catActivityHab2(pass.pass_Dh1_Dh2,:),2)]);
                MeanFiring_120_Dnh1Dh2 = cat(1,MeanFiring_120_Dnh1Dh2,[mean(catActivityHab1(pass.pass_Dnh1_Dh2,:),2),mean(catActivityHab2(pass.pass_Dnh1_Dh2,:),2)]);
                MeanFiring_120_Dh1Dnh2 = cat(1,MeanFiring_120_Dh1Dnh2,[mean(catActivityHab1(pass.pass_Dh1_Dnh2,:),2),mean(catActivityHab2(pass.pass_Dh1_Dnh2,:),2)]);
                MeanFiring_120_Dnh1Dnh2 = cat(1,MeanFiring_120_Dnh1Dnh2,[mean(catActivityHab1(pass.pass_Dnh1_Dnh2,:),2),mean(catActivityHab2(pass.pass_Dnh1_Dnh2,:),2)]);
            end
        end
        
        title([mousename ' Population Mean Interaction Firing Rate across days'])
        saveas(gcf,['PopulationInteractionMeanFiringAcrossTrials'])
        saveas(gcf,['PopulationInteractionMeanFiringAcrossTrials.jpeg'])                              
        
        for i = 1 : length(PassedShuffledCutoff1(:,1))
            if length(PassedShuffledCutoff1(1,:)) > 0
                perInteractions1(i) = length(find(PassedShuffledCutoff1(i,:)))/length(PassedShuffledCutoff1(1,:));
            else
                perInteractions1(i) = 0;
            end
            if length(PassedShuffledCutoff2(1,:)) > 0
                perInteractions2(i) = length(find(PassedShuffledCutoff2(i,:)))/length(PassedShuffledCutoff2(1,:));
            else
                perInteractions2(i) = 0;
            end
            if length(PassedShuffledCutoffHab1(1,:)) > 0
                perInteractionsHab1(i) = length(find(PassedShuffledCutoffHab1(i,:)))/length(PassedShuffledCutoffHab1(1,:));
            else 
                perInteractionsHab1(i) = 0;
            end
            if length(PassedShuffledCutoffHab2(1,:)) > 0
                perInteractionsHab2(i) = length(find(PassedShuffledCutoffHab2(i,:)))/length(PassedShuffledCutoffHab2(1,:));
            else 
                perInteractionsHab2(i) = 0;
            end
            
            if length(PassedShuffledCutoffZ1(1,:)) > 0
                perInteractionsZ1(i) = length(find(PassedShuffledCutoffZ1(i,:)))/length(PassedShuffledCutoffZ1(1,:));
            else
                perInteractionsZ1(i) = 0;
            end
            if length(PassedShuffledCutoffZ2(1,:)) > 0
                perInteractionsZ2(i) = length(find(PassedShuffledCutoffZ2(i,:)))/length(PassedShuffledCutoffZ2(1,:));
            else
                perInteractionsZ2(i) = 0;
            end
            if length(PassedShuffledCutoffZHab1(1,:)) > 0
                perInteractionsZHab1(i) = length(find(PassedShuffledCutoffZHab1(i,:)))/length(PassedShuffledCutoffZHab1(1,:));
            else
                perInteractionsZHab1(i) = 0;
            end
            if length(PassedShuffledCutoffZHab2(1,:)) > 0
                perInteractionsZHab2(i) = length(find(PassedShuffledCutoffZHab2(i,:)))/length(PassedShuffledCutoffZHab2(1,:));
            else
                perInteractionsZHab2(i) = 0;
            end
            
            if length(PassedShuffledCutoffZ1_2(1,:)) > 0
                perInteractionsZ1_2(i) = length(find(PassedShuffledCutoffZ1_2(i,:)))/length(PassedShuffledCutoffZ1_2(1,:));
            else
                perInteractionsZ1_2(i) = 0;
            end
            if length(PassedShuffledCutoffZ2_2(1,:)) > 0
                perInteractionsZ2_2(i) = length(find(PassedShuffledCutoffZ2_2(i,:)))/length(PassedShuffledCutoffZ2_2(1,:));
            else
                perInteractionsZ2_2(i) = 0;
            end
            if length(PassedShuffledCutoffZHab1_2(1,:)) > 0
                perInteractionsZHab1_2(i) = length(find(PassedShuffledCutoffZHab1_2(i,:)))/length(PassedShuffledCutoffZHab1_2(1,:));
            else
                perInteractionsZHab1_2(i) = 0;
            end
            if length(PassedShuffledCutoffZHab2_2(1,:)) > 0
                perInteractionsZHab2_2(i) = length(find(PassedShuffledCutoffZHab2_2(i,:)))/length(PassedShuffledCutoffZHab2_2(1,:));
            else
                perInteractionsZHab2_2(i) = 0;
            end
        end
        for i = 1 : length(PassedShuffledCutoff1prev(:,1))
            if length(PassedShuffledCutoff1prev) > 0
                perInteractions1prev(i) = length(find(PassedShuffledCutoff1prev(i,:)))/length(PassedShuffledCutoff1prev(1,:));
            else 
                perInteractions1prev(i) = 0;
            end
            if length(PassedShuffledCutoff2prev) > 0
                perInteractions2prev(i) = length(find(PassedShuffledCutoff2prev(i,:)))/length(PassedShuffledCutoff2prev(1,:));
            else
                perInteractions2prev(i) = 0;
            end
            if length(PassedShuffledCutoffHab1prev) > 0
                perInteractionsHab1prev(i) = length(find(PassedShuffledCutoffHab1prev(i,:)))/length(PassedShuffledCutoffHab1prev(1,:));
            else
                perInteractionsHab1prev(i) = 0;
            end
            if length(PassedShuffledCutoffHab2prev) > 0
                perInteractionsHab2prev(i) = length(find(PassedShuffledCutoffHab2prev(i,:)))/length(PassedShuffledCutoffHab2prev(1,:));
            else
                perInteractionsHab2prev(i) = 0;
            end
            
            if length(PassedShuffledCutoffZ1prev) > 0
                perInteractionsZ1prev(i) = length(find(PassedShuffledCutoffZ1prev(i,:)))/length(PassedShuffledCutoffZ1prev(1,:));
            else
                perInteractionsZ1prev(i) = 0;
            end
            if length(PassedShuffledCutoffZ2prev) > 0
                perInteractionsZ2prev(i) = length(find(PassedShuffledCutoffZ2prev(i,:)))/length(PassedShuffledCutoffZ2prev(1,:));
            else
                perInteractionsZ2prev(i) = 0;
            end
            if length(PassedShuffledCutoffZHab1prev) > 0
                perInteractionsZHab1prev(i) = length(find(PassedShuffledCutoffZHab1prev(i,:)))/length(PassedShuffledCutoffZHab1prev(1,:));
            else
                perInteractionsZHab1prev(i) = 0;
            end
            if length(PassedShuffledCutoffZHab2prev) > 0
                perInteractionsZHab2prev(i) = length(find(PassedShuffledCutoffZHab2prev(i,:)))/length(PassedShuffledCutoffZHab2prev(1,:));
            else
                perInteractionsZHab2prev(i) = 0;
            end
            
            
            if length(PassedShuffledCutoffZ1prev_2) > 0
                perInteractionsZ1prev_2(i) = length(find(PassedShuffledCutoffZ1prev_2(i,:)))/length(PassedShuffledCutoffZ1prev_2(1,:));
            else
                perInteractionsZ1prev_2(i) = 0;
            end
            if length(PassedShuffledCutoffZ2prev_2) > 0
                perInteractionsZ2prev_2(i) = length(find(PassedShuffledCutoffZ2prev_2(i,:)))/length(PassedShuffledCutoffZ2prev_2(1,:));
            else
                perInteractionsZ2prev_2(i) = 0;
            end
            if length(PassedShuffledCutoffZHab1prev_2) > 0
                perInteractionsZHab1prev_2(i) = length(find(PassedShuffledCutoffZHab1prev_2(i,:)))/length(PassedShuffledCutoffZHab1prev_2(1,:));
            else
                perInteractionsZHab1prev_2(i) = 0;
            end
            if length(PassedShuffledCutoffZHab2prev_2) > 0
                perInteractionsZHab2prev_2(i) = length(find(PassedShuffledCutoffZHab2prev_2(i,:)))/length(PassedShuffledCutoffZHab2prev_2(1,:));
            else
                perInteractionsZHab2prev_2(i) = 0;
            end
        end
        
        if SITnovelty_mins{SITseshnum,2} == 30
            minind = 0;
        else
            minind = 2;
        end
        if contains(SITnovelty_mins{SITseshnum,1},'Different')
            perInteractions{mousecount,1,minind+2} = perInteractions1;
            perInteractions{mousecount,2,minind+2} = perInteractions2;
            perInteractions{mousecount,3,minind+2} = perInteractionsHab1;
            perInteractions{mousecount,4,minind+2} = perInteractionsHab2;
            
            perInteractions{mousecount,1,minind+1} = perInteractions1prev;
            perInteractions{mousecount,2,minind+1} = perInteractions2prev;
            perInteractions{mousecount,3,minind+1} = perInteractionsHab1prev;
            perInteractions{mousecount,4,minind+1} = perInteractionsHab2prev;
            
            perInteractionsZ{mousecount,1,minind+2} = perInteractionsZ1;
            perInteractionsZ{mousecount,2,minind+2} = perInteractionsZ2;
            perInteractionsZ{mousecount,3,minind+2} = perInteractionsZHab1;
            perInteractionsZ{mousecount,4,minind+2} = perInteractionsZHab2;
            
            perInteractionsZ{mousecount,1,minind+1} = perInteractionsZ1prev;
            perInteractionsZ{mousecount,2,minind+1} = perInteractionsZ2prev;
            perInteractionsZ{mousecount,3,minind+1} = perInteractionsZHab1prev;
            perInteractionsZ{mousecount,4,minind+1} = perInteractionsZHab2prev;
                        
            perInteractionsZ_2{mousecount,1,minind+2} = perInteractionsZ1_2;
            perInteractionsZ_2{mousecount,2,minind+2} = perInteractionsZ2_2;
            perInteractionsZ_2{mousecount,3,minind+2} = perInteractionsZHab1_2;
            perInteractionsZ_2{mousecount,4,minind+2} = perInteractionsZHab2_2;
            
            perInteractionsZ_2{mousecount,1,minind+1} = perInteractionsZ1prev_2;
            perInteractionsZ_2{mousecount,2,minind+1} = perInteractionsZ2prev_2;
            perInteractionsZ_2{mousecount,3,minind+1} = perInteractionsZHab1prev_2;
            perInteractionsZ_2{mousecount,4,minind+1} = perInteractionsZHab2prev_2;
        else
            perInteractions{mousecount,1,minind+1} = perInteractions1;
            perInteractions{mousecount,2,minind+1} = perInteractions2;
            perInteractions{mousecount,3,minind+1} = perInteractionsHab1;
            perInteractions{mousecount,4,minind+1} = perInteractionsHab2;
            
            perInteractions{mousecount,1,minind+2} = perInteractions1prev;
            perInteractions{mousecount,2,minind+2} = perInteractions2prev;
            perInteractions{mousecount,3,minind+2} = perInteractionsHab1prev;
            perInteractions{mousecount,4,minind+2} = perInteractionsHab2prev;
            
            perInteractionsZ{mousecount,1,minind+1} = perInteractions1;
            perInteractionsZ{mousecount,2,minind+1} = perInteractions2;
            perInteractionsZ{mousecount,3,minind+1} = perInteractionsHab1;
            perInteractionsZ{mousecount,4,minind+1} = perInteractionsHab2;
            
            perInteractionsZ{mousecount,1,minind+2} = perInteractionsZ1prev;
            perInteractionsZ{mousecount,2,minind+2} = perInteractionsZ2prev;
            perInteractionsZ{mousecount,3,minind+2} = perInteractionsZHab1prev;
            perInteractionsZ{mousecount,4,minind+2} = perInteractionsZHab2prev;
            
            perInteractionsZ_2{mousecount,1,minind+1} = perInteractionsZ1_2;
            perInteractionsZ_2{mousecount,2,minind+1} = perInteractionsZ2_2;
            perInteractionsZ_2{mousecount,3,minind+1} = perInteractionsZHab1_2;
            perInteractionsZ_2{mousecount,4,minind+1} = perInteractionsZHab2_2;
            
            perInteractionsZ_2{mousecount,1,minind+2} = perInteractionsZ1prev_2;
            perInteractionsZ_2{mousecount,2,minind+2} = perInteractionsZ2prev_2;
            perInteractionsZ_2{mousecount,3,minind+2} = perInteractionsZHab1prev_2;
            perInteractionsZ_2{mousecount,4,minind+2} = perInteractionsZHab2prev_2;
        end
        mousecount = mousecount + 1;
        save('InteractionsShufflingPassed.mat','PassedShuffledCutoff1','PassedShuffledCutoff2','PassedShuffledCutoffHab1','PassedShuffledCutoffHab1','PassedShuffledCutoff1prev','PassedShuffledCutoff2prev','PassedShuffledCutoffHab1prev','PassedShuffledCutoffHab1prev',...
            'PassedShuffledCutoffZ1','PassedShuffledCutoffZ2','PassedShuffledCutoffZHab1','PassedShuffledCutoffZHab1','PassedShuffledCutoffZ1prev','PassedShuffledCutoffZ2prev','PassedShuffledCutoffZHab1prev','PassedShuffledCutoffZHab1prev');
        
        
        %% Registration Across days analysis
        load([SITproximityData.folderpaths{1,previousSesh(1)} '\..\..\CellReg\\ms1.mat'],'alignment');
        prevInteractionCal = load([cdold '\' seshdatePrev '\' mousename '\InteractionCalcium.mat']);   
        prevZCal = load([cdold '\' seshdatePrev '\' mousename '\zdataUncut.mat']);   
        if ~isempty(alignment)
            regmap = alignment.alignmentMap{1,2};
            regind = regmap;
            regind(regind>0) = 1;
            regind = regmap(find(sum(regind,2)==2),:);
            Originalregind = [];
            newregind = regind;
            for i = 1 : length(regind(:,1))
                Originalregind(i,2) = gcellind(newregind(i,2));
                Originalregind(i,1) = prevCellinds.gcellind(newregind(i,1));
            end
            meanFR1_1Raw = nansum(prevInteractionCal.catActivity1Raw(regind(:,1),:),2)./length(prevInteractionCal.catActivity1Raw(1,:));
            meanFR2_1Raw = nansum(prevInteractionCal.catActivity2Raw(regind(:,1),:),2)./length(prevInteractionCal.catActivity2Raw(1,:));
            meanFRHab1_1Raw = nansum(prevInteractionCal.catActivityHab1Raw(regind(:,1),:),2)./length(prevInteractionCal.catActivityHab1Raw(1,:));
            meanFRHab2_1Raw = nansum(prevInteractionCal.catActivityHab2Raw(regind(:,1),:),2)./length(prevInteractionCal.catActivityHab2Raw(1,:));
            
            meanFR1_2Raw = nansum(catActivity1Raw(regind(:,2),:),2)./length(catActivity1Raw(1,:));
            meanFR2_2Raw = nansum(catActivity2Raw(regind(:,2),:),2)./length(catActivity2Raw(1,:));
            meanFRHab1_2Raw = nansum(catActivityHab1Raw(regind(:,2),:),2)./length(catActivityHab1Raw(1,:));
            meanFRHab2_2Raw = nansum(catActivityHab2Raw(regind(:,2),:),2)./length(catActivityHab2Raw(1,:));
            
            meanFR1_1z = nansum(prevInteractionCal.catActivity1(regind(:,1),:),2)./length(prevInteractionCal.catActivity1(1,:));
            meanFR2_1z = nansum(prevInteractionCal.catActivity2(regind(:,1),:),2)./length(prevInteractionCal.catActivity2(1,:));
            meanFRHab1_1z = nansum(prevInteractionCal.catActivityHab1(regind(:,1),:),2)./length(prevInteractionCal.catActivityHab1(1,:));
            meanFRHab2_1z = nansum(prevInteractionCal.catActivityHab2(regind(:,1),:),2)./length(prevInteractionCal.catActivityHab2(1,:));
            
            meanFR1_2z = nansum(catActivity1(regind(:,2),:),2)./length(catActivity1(1,:));
            meanFR2_2z = nansum(catActivity2(regind(:,2),:),2)./length(catActivity2(1,:));
            meanFRHab1_2z = nansum(catActivityHab1(regind(:,2),:),2)./length(catActivityHab1(1,:));
            meanFRHab2_2z = nansum(catActivityHab2(regind(:,2),:),2)./length(catActivityHab2(1,:));
            
            maxFRraw = max(cat(1,meanFR1_1Raw,meanFR1_2Raw,meanFR2_1Raw,meanFR2_2Raw,meanFRHab1_1Raw,meanFRHab1_2Raw,meanFRHab2_1Raw,meanFRHab2_2Raw));
            minFRraw = min(cat(1,meanFR1_1Raw,meanFR1_2Raw,meanFR2_1Raw,meanFR2_2Raw,meanFRHab1_1Raw,meanFRHab1_2Raw,meanFRHab2_1Raw,meanFRHab2_2Raw));
            
            maxFRz = max(cat(1,meanFR1_1z,meanFR1_2z,meanFR2_1z,meanFR2_2z,meanFRHab1_1z,meanFRHab1_2z,meanFRHab2_1z,meanFRHab2_2z));
            minFRz = min(cat(1,meanFR1_1z,meanFR1_2z,meanFR2_1z,meanFR2_2z,meanFRHab1_1z,meanFRHab1_2z,meanFRHab2_1z,meanFRHab2_2z));
                        
            [~, ind2]= intersect(regind(:,2),unique(cat(1,find(p1z),find(p2z),find(ph1z),find(ph2z))));
            [~, ind1] = intersect(regind(:,1),unique(cat(1,find(p1prevz),find(p2prevz),find(ph1prevz),find(ph2prevz))));
            
            ind3 = unique(cat(1,ind2,ind1));
            regindpassed = regind(ind3,:);
            
            meanFR1_1z_passed_CR = nansum(prevInteractionCal.catActivity1(regindpassed(:,1),:),2)./length(prevInteractionCal.catActivity1(1,:));
            meanFR2_1z_passed_CR = nansum(prevInteractionCal.catActivity2(regindpassed(:,1),:),2)./length(prevInteractionCal.catActivity2(1,:));
            meanFRHab1_1z_passed_CR = nansum(prevInteractionCal.catActivityHab1(regindpassed(:,1),:),2)./length(prevInteractionCal.catActivityHab1(1,:));
            meanFRHab2_1z_passed_CR = nansum(prevInteractionCal.catActivityHab2(regindpassed(:,1),:),2)./length(prevInteractionCal.catActivityHab2(1,:));
            
            meanFR1_2z_passed_CR = nansum(catActivity1(regindpassed(:,2),:),2)./length(catActivity1(1,:));
            meanFR2_2z_passed_CR = nansum(catActivity2(regindpassed(:,2),:),2)./length(catActivity2(1,:));
            meanFRHab1_2z_passed_CR = nansum(catActivityHab1(regindpassed(:,2),:),2)./length(catActivityHab1(1,:));
            meanFRHab2_2z_passed_CR = nansum(catActivityHab2(regindpassed(:,2),:),2)./length(catActivityHab2(1,:));
            
            maxFRz_passed = max(cat(1,meanFR1_1z_passed_CR,meanFR1_2z_passed_CR,meanFR2_1z_passed_CR,meanFR2_2z_passed_CR,meanFRHab1_1z_passed_CR,meanFRHab1_2z_passed_CR,meanFRHab2_1z_passed_CR,meanFRHab2_2z_passed_CR));
            minFRz_passed = min(cat(1,meanFR1_1z_passed_CR,meanFR1_2z_passed_CR,meanFR2_1z_passed_CR,meanFR2_2z_passed_CR,meanFRHab1_1z_passed_CR,meanFRHab1_2z_passed_CR,meanFRHab2_1z_passed_CR,meanFRHab2_2z_passed_CR));
            
            % non Cell Reg 
            meanFR1_1RawAll = nansum(prevInteractionCal.catActivity1Raw,2)./length(prevInteractionCal.catActivity1Raw(1,:));
            meanFR2_1RawAll = nansum(prevInteractionCal.catActivity2Raw,2)./length(prevInteractionCal.catActivity2Raw(1,:));
            meanFRHab1_1RawAll = nansum(prevInteractionCal.catActivityHab1Raw,2)./length(prevInteractionCal.catActivityHab1Raw(1,:));
            meanFRHab2_1RawAll = nansum(prevInteractionCal.catActivityHab2Raw,2)./length(prevInteractionCal.catActivityHab2Raw(1,:));
            
            meanFR1_2RawAll = nansum(catActivity1Raw,2)./length(catActivity1Raw(1,:));
            meanFR2_2RawAll = nansum(catActivity2Raw,2)./length(catActivity2Raw(1,:));
            meanFRHab1_2RawAll = nansum(catActivityHab1Raw,2)./length(catActivityHab1Raw(1,:));
            meanFRHab2_2RawAll = nansum(catActivityHab2Raw,2)./length(catActivityHab2Raw(1,:));
            
            meanFR1_1zAll = nansum(prevInteractionCal.catActivity1,2)./length(prevInteractionCal.catActivity1(1,:));
            meanFR2_1zAll = nansum(prevInteractionCal.catActivity2,2)./length(prevInteractionCal.catActivity2(1,:));
            meanFRHab1_1zAll = nansum(prevInteractionCal.catActivityHab1,2)./length(prevInteractionCal.catActivityHab1(1,:));
            meanFRHab2_1zAll = nansum(prevInteractionCal.catActivityHab2,2)./length(prevInteractionCal.catActivityHab2(1,:));
            
            meanFR1_2zAll = nansum(catActivity1,2)./length(catActivity1(1,:));
            meanFR2_2zAll = nansum(catActivity2,2)./length(catActivity2(1,:));
            meanFRHab1_2zAll = nansum(catActivityHab1,2)./length(catActivityHab1(1,:));
            meanFRHab2_2zAll = nansum(catActivityHab2,2)./length(catActivityHab2(1,:));
            
            maxFRrawAll = max(cat(1,meanFR1_1RawAll,meanFR1_2RawAll,meanFR2_1RawAll,meanFR2_2RawAll,meanFRHab1_1RawAll,meanFRHab1_2RawAll,meanFRHab2_1RawAll,meanFRHab2_2RawAll));
            minFRrawAll = min(cat(1,meanFR1_1RawAll,meanFR1_2RawAll,meanFR2_1RawAll,meanFR2_2RawAll,meanFRHab1_1RawAll,meanFRHab1_2RawAll,meanFRHab2_1RawAll,meanFRHab2_2RawAll));
            
            maxFRzAll = max(cat(1,meanFR1_1zAll,meanFR1_2zAll,meanFR2_1zAll,meanFR2_2zAll,meanFRHab1_1zAll,meanFRHab1_2zAll,meanFRHab2_1zAll,meanFRHab2_2zAll));
            minFRzAll = min(cat(1,meanFR1_1zAll,meanFR1_2zAll,meanFR2_1zAll,meanFR2_2zAll,meanFRHab1_1zAll,meanFRHab1_2zAll,meanFRHab2_1zAll,meanFRHab2_2zAll));            
            
            meanFR1_1zAll_passed = nansum(prevInteractionCal.catActivity1(p1prevz,:),2)./length(prevInteractionCal.catActivity1(1,:));
            meanFR2_1zAll_passed = nansum(prevInteractionCal.catActivity2(p2prevz,:),2)./length(prevInteractionCal.catActivity2(1,:));
            meanFRHab1_1zAll_passed = nansum(prevInteractionCal.catActivityHab1(ph1prevz,:),2)./length(prevInteractionCal.catActivityHab1(1,:));
            meanFRHab2_1zAll_passed = nansum(prevInteractionCal.catActivityHab2(ph2prevz,:),2)./length(prevInteractionCal.catActivityHab2(1,:));
            
            meanFR1_2zAll_passed = nansum(catActivity1(p1z,:),2)./length(catActivity1(1,:));
            meanFR2_2zAll_passed = nansum(catActivity2(p2z,:),2)./length(catActivity2(1,:));
            meanFRHab1_2zAll_passed = nansum(catActivityHab1(ph1z,:),2)./length(catActivityHab1(1,:));
            meanFRHab2_2zAll_passed = nansum(catActivityHab2(ph2z,:),2)./length(catActivityHab2(1,:));                        
            
            maxFRzAll_passed = max(cat(1,meanFR1_1zAll_passed,meanFR1_2zAll_passed,meanFR2_1zAll_passed,meanFR2_2zAll_passed,meanFRHab1_1zAll_passed,meanFRHab1_2zAll_passed,meanFRHab2_1zAll_passed,meanFRHab2_2zAll_passed));
            minFRzAll_passed = min(cat(1,meanFR1_1zAll_passed,meanFR1_2zAll_passed,meanFR2_1zAll_passed,meanFR2_2zAll_passed,meanFRHab1_1zAll_passed,meanFRHab1_2zAll_passed,meanFRHab2_1zAll_passed,meanFRHab2_2zAll_passed));
            
            temp1 = zeros(length(catActivity1(:,1)),1);
            temp2 = zeros(length(catActivity2(:,1)),1);
            temp3 = zeros(length(catActivityHab1(:,1)),1);
            temp4 = zeros(length(catActivityHab2(:,1)),1);

            p1z_CR = intersect(regind(:,2),p1z_CR);
            p2z_CR = intersect(regind(:,2),p2z_CR);
            ph1z_CR = intersect(regind(:,2),ph1z_CR);
            ph2z_CR = intersect(regind(:,2),ph2z_CR);
            
            p1prevz_CR = intersect(regind(:,1),p1prevz_CR);
            p2prevz_CR = intersect(regind(:,1),p2prevz_CR);
            ph1prevz_CR = intersect(regind(:,1),ph1prevz_CR);
            ph2prevz_CR = intersect(regind(:,1),ph2prevz_CR);
            
            [~,ind_1reg] = intersect(regind(:,2),p1z_CR);
            [~,ind_2reg] = intersect(regind(:,2),p2z_CR);
            [~,ind_h1reg] = intersect(regind(:,2),ph1z_CR);
            [~,ind_h2reg] = intersect(regind(:,2),ph2z_CR);
            [~,ind_1regprev] = intersect(regind(:,1),p1prevz_CR);
            [~,ind_2regprev] = intersect(regind(:,1),p2prevz_CR);
            [~,ind_h1regprev] = intersect(regind(:,1),ph1prevz_CR);
            [~,ind_h2regprev] = intersect(regind(:,1),ph2prevz_CR);
            
            [~,ind_1reg_not,~] = setxor(regind(:,2),p1z_CR);
            [~,ind_2reg_not,~] = setxor(regind(:,2),p2z_CR);
            [~,ind_h1reg_not,~] = setxor(regind(:,2),ph1z_CR);
            [~,ind_h2reg_not,~] = setxor(regind(:,2),ph2z_CR);
            [~,ind_1regprev_not,~] = setxor(regind(:,1),p1prevz_CR);
            [~,ind_2regprev_not,~] = setxor(regind(:,1),p2prevz_CR);
            [~,ind_h1regprev_not,~] = setxor(regind(:,1),ph1prevz_CR);
            [~,ind_h2regprev_not,~] = setxor(regind(:,1),ph2prevz_CR);    
            
            pass = [];
            
            if SITnovelty_mins{SITseshnum,2} == 30
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    pass = passCellReg(regind,ind_1reg,ind_2reg,ind_h1reg,ind_h2reg,ind_1reg_not,ind_2reg_not,ind_h1reg_not,ind_h2reg_not,ind_1regprev,ind_2regprev,ind_h1regprev,ind_h2regprev,ind_1regprev_not,ind_2regprev_not,ind_h1regprev_not,ind_h2regprev_not);
                    %Trial1 Same vs Trial1 Diff
                    MeanFiringCellReg_30_S1D1 = cat(1,MeanFiringCellReg_30_S1D1,[mean(catActivity1(pass.pass_S1_D1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_S1_D1(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn1D1 = cat(1,MeanFiringCellReg_30_Sn1D1,[mean(catActivity1(pass.pass_Sn1_D1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_Sn1_D1(:,1),:),2)]);
                    MeanFiringCellReg_30_S1Dn1 = cat(1,MeanFiringCellReg_30_S1Dn1,[mean(catActivity1(pass.pass_S1_Dn1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_S1_Dn1(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn1Dn1 = cat(1,MeanFiringCellReg_30_Sn1Dn1,[mean(catActivity1(pass.pass_Sn1_Dn1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_Sn1_Dn1(:,1),:),2)]);
                    %Trial2 Same vs Trial1 Diff
                    MeanFiringCellReg_30_S2D1 = cat(1,MeanFiringCellReg_30_S2D1,[mean(catActivity2(pass.pass_S2_D1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_S2_D1(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn2D1 = cat(1,MeanFiringCellReg_30_Sn2D1,[mean(catActivity2(pass.pass_Sn2_D1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_Sn2_D1(:,1),:),2)]);
                    MeanFiringCellReg_30_S2Dn1 = cat(1,MeanFiringCellReg_30_S2Dn1,[mean(catActivity2(pass.pass_S2_Dn1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_S2_Dn1(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn2Dn1 = cat(1,MeanFiringCellReg_30_Sn2Dn1,[mean(catActivity2(pass.pass_Sn2_Dn1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_Sn2_Dn1(:,1),:),2)]);
                    %Trial2 Same vs Trial2 Diff
                    MeanFiringCellReg_30_S2D2 = cat(1,MeanFiringCellReg_30_S2D2,[mean(catActivity2(pass.pass_S2_D2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_S2_D2(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn2D2 = cat(1,MeanFiringCellReg_30_Sn2D2,[mean(catActivity2(pass.pass_Sn2_D2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_Sn2_D2(:,1),:),2)]);
                    MeanFiringCellReg_30_S2Dn2 = cat(1,MeanFiringCellReg_30_S2Dn2,[mean(catActivity2(pass.pass_S2_Dn2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_S2_Dn2(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn2Dn2 = cat(1,MeanFiringCellReg_30_Sn2Dn2,[mean(catActivity2(pass.pass_Sn2_Dn2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_Sn2_Dn2(:,1),:),2)]);
                    %Trial1 Same vs Trial2 Diff
                    MeanFiringCellReg_30_S1D2 = cat(1,MeanFiringCellReg_30_S1D2,[mean(catActivity1(pass.pass_S1_D2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_S1_D2(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn1D2 = cat(1,MeanFiringCellReg_30_Sn1D2,[mean(catActivity1(pass.pass_Sn1_D2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_Sn1_D2(:,1),:),2)]);
                    MeanFiringCellReg_30_S1Dn2 = cat(1,MeanFiringCellReg_30_S1Dn2,[mean(catActivity1(pass.pass_S1_Dn2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_S1_Dn2(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn1Dn2 = cat(1,MeanFiringCellReg_30_Sn1Dn2,[mean(catActivity1(pass.pass_Sn1_Dn2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_Sn1_Dn2(:,1),:),2)]);
                    %Hab1 Same vs Hab1 Diff
                    MeanFiringCellReg_30_Sh1Dh1 = cat(1,MeanFiringCellReg_30_Sh1Dh1,[mean(catActivityHab1(pass.pass_Sh1_Dh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Sh1_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_30_Snh1Dh1 = cat(1,MeanFiringCellReg_30_Snh1Dh1,[mean(catActivityHab1(pass.pass_Sh1_Dnh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Sh1_Dnh1(:,1),:),2)]);
                    MeanFiringCellReg_30_Sh1Dnh1 = cat(1,MeanFiringCellReg_30_Sh1Dnh1,[mean(catActivityHab1(pass.pass_Snh1_Dh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Snh1_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_30_Snh1Dnh1 = cat(1,MeanFiringCellReg_30_Snh1Dnh1,[mean(catActivityHab1(pass.pass_Snh1_Dnh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Snh1_Dnh1(:,1),:),2)]);
                    %Hab1 Same vs Hab2 Diff
                    MeanFiringCellReg_30_Sh1Dh2 = cat(1,MeanFiringCellReg_30_Sh1Dh2,[mean(catActivityHab1(pass.pass_Sh1_Dh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Sh1_Dh2(:,1),:),2)]);
                    MeanFiringCellReg_30_Snh1Dh2 = cat(1,MeanFiringCellReg_30_Snh1Dh2,[mean(catActivityHab1(pass.pass_Snh1_Dh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Snh1_Dh2(:,1),:),2)]);
                    MeanFiringCellReg_30_Sh1Dnh2 = cat(1,MeanFiringCellReg_30_Sh1Dnh2,[mean(catActivityHab1(pass.pass_Sh1_Dnh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Sh1_Dnh2(:,1),:),2)]);
                    MeanFiringCellReg_30_Snh1Dnh2 = cat(1,MeanFiringCellReg_30_Snh1Dnh2,[mean(catActivityHab1(pass.pass_Snh1_Dnh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Snh1_Dnh2(:,1),:),2)]);
                    %Hab2 Same vs Hab2 Diff
                    MeanFiringCellReg_30_Sh2Dh2 = cat(1,MeanFiringCellReg_30_Sh2Dh2,[mean(catActivityHab2(pass.pass_Sh2_Dh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Sh2_Dh2(:,1),:),2)]);
                    MeanFiringCellReg_30_Snh2Dh2 = cat(1,MeanFiringCellReg_30_Snh2Dh2,[mean(catActivityHab2(pass.pass_Snh2_Dh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Snh2_Dh2(:,1),:),2)]);
                    MeanFiringCellReg_30_Sh2Dnh2 = cat(1,MeanFiringCellReg_30_Sh2Dnh2,[mean(catActivityHab2(pass.pass_Sh2_Dnh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Sh2_Dnh2(:,1),:),2)]);
                    MeanFiringCellReg_30_Snh2Dnh2 = cat(1,MeanFiringCellReg_30_Snh2Dnh2,[mean(catActivityHab2(pass.pass_Snh2_Dnh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Snh2_Dnh2(:,1),:),2)]);
                    %Hab2 Same vs Hab1 Diff
                    MeanFiringCellReg_30_Sh2Dh1 = cat(1,MeanFiringCellReg_30_Sh2Dh1,[mean(catActivityHab2(pass.pass_Sh2_Dh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Sh2_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_30_Snh2Dh1 = cat(1,MeanFiringCellReg_30_Snh2Dh1,[mean(catActivityHab2(pass.pass_Snh2_Dh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Snh2_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_30_Sh2Dnh1 = cat(1,MeanFiringCellReg_30_Sh2Dnh1,[mean(catActivityHab2(pass.pass_Sh2_Dnh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Sh2_Dnh1(:,1),:),2)]);
                    MeanFiringCellReg_30_Snh2Dnh1 = cat(1,MeanFiringCellReg_30_Snh2Dnh1,[mean(catActivityHab2(pass.pass_Snh2_Dnh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Snh2_Dnh1(:,1),:),2)]);
                    %Trial1 Same vs Trial2 Same
                    MeanFiringCellReg_30_S1S2 = cat(1,MeanFiringCellReg_30_S1D1,[mean(catActivity1(pass.pass_S1_S2(:,2),:),2),mean(catActivity2(pass.pass_S1_S2(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn1S2 = cat(1,MeanFiringCellReg_30_Sn1D1,[mean(catActivity1(pass.pass_Sn1_S2(:,2),:),2),mean(catActivity2(pass.pass_Sn1_S2(:,2),:),2)]);
                    MeanFiringCellReg_30_S1Sn2 = cat(1,MeanFiringCellReg_30_S1Dn1,[mean(catActivity1(pass.pass_S1_Sn2(:,2),:),2),mean(catActivity2(pass.pass_S1_Sn2(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn1Sn2 = cat(1,MeanFiringCellReg_30_Sn1Dn1,[mean(catActivity1(pass.pass_Sn1_Sn2(:,2),:),2),mean(catActivity2(pass.pass_Sn1_Sn2(:,2),:),2)]);
                    %Trial1 Diff vs Trial2 Diff
                    MeanFiringCellReg_30_D1D2 = cat(1,MeanFiringCellReg_30_D1D2,[mean(prev4b.catActivity1(pass.pass_D1_D2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_D1_D2(:,1),:),2)]);
                    MeanFiringCellReg_30_Dn1D2 = cat(1,MeanFiringCellReg_30_Dn1D2,[mean(prev4b.catActivity1(pass.pass_Dn1_D2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_Dn1_D2(:,1),:),2)]);
                    MeanFiringCellReg_30_D1Dn2 = cat(1,MeanFiringCellReg_30_D1Dn2,[mean(prev4b.catActivity1(pass.pass_D1_Dn2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_D1_Dn2(:,1),:),2)]);
                    MeanFiringCellReg_30_Dn1Dn2 = cat(1,MeanFiringCellReg_30_Dn1Dn2,[mean(prev4b.catActivity1(pass.pass_Dn1_Dn2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_Dn1_Dn2(:,1),:),2)]);
                    
                    
                    if isempty(Chrono30Hab1Same)
                        Chrono30Hab1Same = chronoall{previousSesh(2),1}(regind(:,2),:);
                        Chrono301Same = chronoall{previousSesh(2),2}(regind(:,2),:);
                        Chrono30Hab2Same = chronoall{previousSesh(2),3}(regind(:,2),:);
                        Chrono302Same = chronoall{previousSesh(2),4}(regind(:,2),:);
                        
                        Chrono30Hab1Diff = chronoall{previousSesh(1),1}(regind(:,1),:);
                        Chrono301Diff = chronoall{previousSesh(1),2}(regind(:,1),:);
                        Chrono30Hab2Diff = chronoall{previousSesh(1),3}(regind(:,1),:);
                        Chrono302Diff = chronoall{previousSesh(1),4}(regind(:,1),:);                        
                    else
                        Chrono30Hab1Same = cat(1,Chrono30Hab1Same,chronoall{previousSesh(2),1}(regind(:,2),:));
                        Chrono301Same = cat(1,Chrono301Same,chronoall{previousSesh(2),2}(regind(:,2),:));
                        Chrono30Hab2Same = cat(1,Chrono30Hab2Same,chronoall{previousSesh(2),3}(regind(:,2),:));
                        Chrono302Same = cat(1,Chrono302Same,chronoall{previousSesh(2),4}(regind(:,2),:));
                        
                        Chrono30Hab1Diff = cat(1,Chrono30Hab1Diff,chronoall{previousSesh(1),1}(regind(:,1),:));
                        Chrono301Diff = cat(1,Chrono301Diff,chronoall{previousSesh(1),2}(regind(:,1),:));
                        Chrono30Hab2Diff = cat(1,Chrono30Hab2Diff,chronoall{previousSesh(1),3}(regind(:,1),:));
                        Chrono302Diff = cat(1,Chrono302Diff,chronoall{previousSesh(1),4}(regind(:,1),:));
                    end                                        
                    % trial 1 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,1})
                        if i > 1 && enterzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i-1) >= secframes2 && enterzonesall{SITseshnum,1}(i) + secframes2 < length(ztrial1(1,:)) && exitzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i) >= secframes2
                            bout1_2sec_30Same_CellReg = cat(1,bout1_2sec_30Same_CellReg,ztrial1(gcellind(regind(:,2)),enterzonesall{SITseshnum,1}(i)-secframes2:enterzonesall{SITseshnum,1}(i)+secframes2));
                        elseif i == 1 && enterzonesall{SITseshnum,1}(1) > secframes2 && exitzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i) >= secframes2
                            bout1_2sec_30Same_CellReg = cat(1,bout1_2sec_30Same_CellReg,ztrial1(gcellind(regind(:,2)),enterzonesall{SITseshnum,1}(i)-secframes2:enterzonesall{SITseshnum,1}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),1})
                        if i > 1 && enterzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i-1) >= secframes2 && enterzonesall{previousSesh(1),1}(i) + secframes2 < length(prevZCal.ztrial1(1,:)) && exitzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i) >= secframes2
                            bout1_2sec_30Diff_CellReg = cat(1,bout1_2sec_30Diff_CellReg,prevZCal.ztrial1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),1}(i)-secframes2:enterzonesall{previousSesh(1),1}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),1}(1) > secframes2 && exitzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i) >= secframes2
                            bout1_2sec_30Diff_CellReg = cat(1,bout1_2sec_30Diff_CellReg,prevZCal.ztrial1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),1}(i)-secframes2:enterzonesall{previousSesh(1),1}(i)+secframes2));
                        end
                    end
                    % trial 2 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,2})
                        if i > 1 && enterzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i-1) >= secframes2 && enterzonesall{SITseshnum,2}(i) + secframes2 < length(ztrial2(1,:)) && exitzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i) >= secframes2
                            bout2_2sec_30Same_CellReg = cat(1,bout2_2sec_30Same_CellReg,ztrial2(gcellind(regind(:,2)),enterzonesall{SITseshnum,2}(i)-secframes2:enterzonesall{SITseshnum,2}(i)+secframes2));
                        elseif i == 1 && enterzonesall{SITseshnum,1}(1) > secframes2 && exitzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i) >= secframes2
                            bout2_2sec_30Same_CellReg = cat(1,bout2_2sec_30Same_CellReg,ztrial2(gcellind(regind(:,2)),enterzonesall{SITseshnum,2}(i)-secframes2:enterzonesall{SITseshnum,2}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),2})
                        if i > 1 && enterzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i-1) >= secframes2 && enterzonesall{previousSesh(1),2}(i) + secframes2 < length(prevZCal.ztrial2(1,:)) && exitzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i) >= secframes2
                            bout2_2sec_30Diff_CellReg = cat(1,bout2_2sec_30Diff_CellReg,prevZCal.ztrial2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),2}(i)-secframes2:enterzonesall{previousSesh(1),2}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),2}(1) > secframes2 && exitzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i) >= secframes2
                            bout2_2sec_30Diff_CellReg = cat(1,bout2_2sec_30Diff_CellReg,prevZCal.ztrial2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),2}(i)-secframes2:enterzonesall{previousSesh(1),2}(i)+secframes2));
                        end
                    end
                    % Hab 1 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,3})
                        if i > 1 && enterzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i-1) >= secframes2 && enterzonesall{SITseshnum,3}(i) + secframes2 < length(zhab1(1,:)) && exitzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i) >= secframes2 
                            boutHab1_2sec_30Same_CellReg = cat(1,boutHab1_2sec_30Same_CellReg,zhab1(gcellind(regind(:,2)),enterzonesall{SITseshnum,3}(i)-secframes2:enterzonesall{SITseshnum,3}(i)+secframes2));
                        elseif i == 1 && enterzonesall{SITseshnum,1}(1) > secframes2 && exitzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i) >= secframes2
                            boutHab1_2sec_30Same_CellReg = cat(1,boutHab1_2sec_30Same_CellReg,zhab1(gcellind(regind(:,2)),enterzonesall{SITseshnum,3}(i)-secframes2:enterzonesall{SITseshnum,3}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),3})
                        if i > 1 && enterzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i-1) >= secframes2 && enterzonesall{previousSesh(1),3}(i) + secframes2 < length(prevZCal.zhab1(1,:)) && exitzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i) >= secframes2
                            boutHab1_2sec_30Diff_CellReg = cat(1,boutHab1_2sec_30Diff_CellReg,prevZCal.zhab1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),3}(i)-secframes2:enterzonesall{previousSesh(1),3}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),3}(1) > secframes2 && exitzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i) >= secframes2
                            boutHab1_2sec_30Diff_CellReg = cat(1,boutHab1_2sec_30Diff_CellReg,prevZCal.zhab1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),3}(i)-secframes2:enterzonesall{previousSesh(1),3}(i)+secframes2));
                        end
                    end
                    % Hab 2 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,4})
                        if i > 1 && enterzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i-1) >= secframes2 && enterzonesall{SITseshnum,4}(i) + secframes2 < length(zhab2(1,:)) && exitzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i) >= secframes2
                            boutHab2_2sec_30Same_CellReg = cat(1,boutHab2_2sec_30Same_CellReg,zhab2(gcellind(regind(:,2)),enterzonesall{SITseshnum,4}(i)-secframes2:enterzonesall{SITseshnum,4}(i)+secframes2));
                        elseif i == 1 && enterzonesall{SITseshnum,1}(1) > secframes2 && exitzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i) >= secframes2
                            boutHab2_2sec_30Same_CellReg = cat(1,boutHab2_2sec_30Same_CellReg,zhab2(gcellind(regind(:,2)),enterzonesall{SITseshnum,4}(i)-secframes2:enterzonesall{SITseshnum,4}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),4})
                        if i > 1 && enterzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i-1) >= secframes2 && enterzonesall{previousSesh(1),4}(i) + secframes2 < length(prevZCal.zhab2(1,:)) && exitzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i) >= secframes2
                            boutHab2_2sec_30Diff_CellReg = cat(1,boutHab2_2sec_30Diff_CellReg,prevZCal.zhab2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),4}(i)-secframes2:enterzonesall{previousSesh(1),4}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),4}(1) > secframes2 && exitzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i) >= secframes2
                            boutHab2_2sec_30Diff_CellReg = cat(1,boutHab2_2sec_30Diff_CellReg,prevZCal.zhab2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),4}(i)-secframes2:enterzonesall{previousSesh(1),4}(i)+secframes2));
                        end
                    end
                else
                    pass = passCellReg(regind,ind_1regprev,ind_2regprev,ind_h1regprev,ind_h2regprev,ind_1regprev_not,ind_2regprev_not,ind_h1regprev_not,ind_h2regprev_not,ind_1reg,ind_2reg,ind_h1reg,ind_h2reg,ind_1reg_not,ind_2reg_not,ind_h1reg_not,ind_h2reg_not);
                    %Trial1 Same vs Trial1 Diff
                    MeanFiringCellReg_30_S1D1 = cat(1,MeanFiringCellReg_30_S1D1,[mean(prev4b.catActivity1(pass.pass_S1_D1(:,1),:),2),mean(catActivity1(pass.pass_S1_D1(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn1D1 = cat(1,MeanFiringCellReg_30_Sn1D1,[mean(prev4b.catActivity1(pass.pass_Sn1_D1(:,1),:),2),mean(catActivity1(pass.pass_Sn1_D1(:,2),:),2)]);
                    MeanFiringCellReg_30_S1Dn1 = cat(1,MeanFiringCellReg_30_S1Dn1,[mean(prev4b.catActivity1(pass.pass_S1_Dn1(:,1),:),2),mean(catActivity1(pass.pass_S1_Dn1(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn1Dn1 = cat(1,MeanFiringCellReg_30_Sn1Dn1,[mean(prev4b.catActivity1(pass.pass_Sn1_Dn1(:,1),:),2),mean(catActivity1(pass.pass_Sn1_Dn1(:,2),:),2)]);
                    %Trial2 Same vs Trial1 Diff
                    MeanFiringCellReg_30_S2D1 = cat(1,MeanFiringCellReg_30_S2D1,[mean(prev4b.catActivity2(pass.pass_S2_D1(:,1),:),2),mean(catActivity1(pass.pass_S2_D1(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn2D1 = cat(1,MeanFiringCellReg_30_Sn2D1,[mean(prev4b.catActivity2(pass.pass_Sn2_D1(:,1),:),2),mean(catActivity1(pass.pass_Sn2_D1(:,2),:),2)]);
                    MeanFiringCellReg_30_S2Dn1 = cat(1,MeanFiringCellReg_30_S2Dn1,[mean(prev4b.catActivity2(pass.pass_S2_Dn1(:,1),:),2),mean(catActivity1(pass.pass_S2_Dn1(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn2Dn1 = cat(1,MeanFiringCellReg_30_Sn2Dn1,[mean(prev4b.catActivity2(pass.pass_Sn2_Dn1(:,1),:),2),mean(catActivity1(pass.pass_Sn2_Dn1(:,2),:),2)]);
                    %Trial2 Same vs Trial2 Diff
                    MeanFiringCellReg_30_S2D2 = cat(1,MeanFiringCellReg_30_S2D2,[mean(prev4b.catActivity2(pass.pass_S2_D2(:,1),:),2),mean(catActivity2(pass.pass_S2_D2(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn2D2 = cat(1,MeanFiringCellReg_30_Sn2D2,[mean(prev4b.catActivity2(pass.pass_Sn2_D2(:,1),:),2),mean(catActivity2(pass.pass_Sn2_D2(:,2),:),2)]);
                    MeanFiringCellReg_30_S2Dn2 = cat(1,MeanFiringCellReg_30_S2Dn2,[mean(prev4b.catActivity2(pass.pass_S2_Dn2(:,1),:),2),mean(catActivity2(pass.pass_S2_Dn2(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn2Dn2 = cat(1,MeanFiringCellReg_30_Sn2Dn2,[mean(prev4b.catActivity2(pass.pass_Sn2_Dn2(:,1),:),2),mean(catActivity2(pass.pass_Sn2_Dn2(:,2),:),2)]);
                    %Trial1 Same vs Trial2 Diff
                    MeanFiringCellReg_30_S1D2 = cat(1,MeanFiringCellReg_30_S1D2,[mean(prev4b.catActivity1(pass.pass_S1_D2(:,1),:),2),mean(catActivity2(pass.pass_S1_D2(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn1D2 = cat(1,MeanFiringCellReg_30_Sn1D2,[mean(prev4b.catActivity1(pass.pass_Sn1_D2(:,1),:),2),mean(catActivity2(pass.pass_Sn1_D2(:,2),:),2)]);
                    MeanFiringCellReg_30_S1Dn2 = cat(1,MeanFiringCellReg_30_S1Dn2,[mean(prev4b.catActivity1(pass.pass_S1_Dn2(:,1),:),2),mean(catActivity2(pass.pass_S1_Dn2(:,2),:),2)]);
                    MeanFiringCellReg_30_Sn1Dn2 = cat(1,MeanFiringCellReg_30_Sn1Dn2,[mean(prev4b.catActivity1(pass.pass_Sn1_Dn2(:,1),:),2),mean(catActivity2(pass.pass_Sn1_Dn2(:,2),:),2)]);
                    %Hab1 Same vs Hab1 Diff
                    MeanFiringCellReg_30_Sh1Dh1 = cat(1,MeanFiringCellReg_30_Sh1Dh1,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Dh1(:,1),:),2),mean(catActivityHab1(pass.pass_Sh1_Dh1(:,2),:),2)]);
                    MeanFiringCellReg_30_Snh1Dh1 = cat(1,MeanFiringCellReg_30_Snh1Dh1,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Dh1(:,1),:),2),mean(catActivityHab1(pass.pass_Snh1_Dh1(:,2),:),2)]);
                    MeanFiringCellReg_30_Sh1Dnh1 = cat(1,MeanFiringCellReg_30_Sh1Dnh1,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Dnh1(:,1),:),2),mean(catActivityHab1(pass.pass_Sh1_Dnh1(:,2),:),2)]);
                    MeanFiringCellReg_30_Snh1Dnh1 = cat(1,MeanFiringCellReg_30_Snh1Dnh1,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Dnh1(:,1),:),2),mean(catActivityHab1(pass.pass_Snh1_Dnh1(:,2),:),2)]);
                    %Hab1 Same vs Hab2 Diff
                    MeanFiringCellReg_30_Sh1Dh2 = cat(1,MeanFiringCellReg_30_Sh1Dh2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Dh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sh1_Dh2(:,2),:),2)]);
                    MeanFiringCellReg_30_Snh1Dh2 = cat(1,MeanFiringCellReg_30_Snh1Dh2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Dh2(:,1),:),2),mean(catActivityHab2(pass.pass_Snh1_Dh2(:,2),:),2)]);
                    MeanFiringCellReg_30_Sh1Dnh2 = cat(1,MeanFiringCellReg_30_Sh1Dnh2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Dnh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sh1_Dnh2(:,2),:),2)]);
                    MeanFiringCellReg_30_Snh1Dnh2 = cat(1,MeanFiringCellReg_30_Snh1Dnh2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Dnh2(:,1),:),2),mean(catActivityHab2(pass.pass_Snh1_Dnh2(:,2),:),2)]);
                    %Hab2 Same vs Hab2 Diff
                    MeanFiringCellReg_30_Sh2Dh2 = cat(1,MeanFiringCellReg_30_Sh2Dh2,[mean(prev4b.catActivityHab2(pass.pass_Sh2_Dh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sh2_Dh2(:,2),:),2)]);
                    MeanFiringCellReg_30_Snh2Dh2 = cat(1,MeanFiringCellReg_30_Snh2Dh2,[mean(prev4b.catActivityHab2(pass.pass_Snh2_Dh2(:,1),:),2),mean(catActivityHab2(pass.pass_Snh2_Dh2(:,2),:),2)]);
                    MeanFiringCellReg_30_Sh2Dnh2 = cat(1,MeanFiringCellReg_30_Sh2Dnh2,[mean(prev4b.catActivityHab2(pass.pass_Sh2_Dnh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sh2_Dnh2(:,2),:),2)]);
                    MeanFiringCellReg_30_Snh2Dnh2 = cat(1,MeanFiringCellReg_30_Snh2Dnh2,[mean(prev4b.catActivityHab2(pass.pass_Snh2_Dnh2(:,1),:),2),mean(catActivityHab2(pass.pass_Snh2_Dnh2(:,2),:),2)]);
                    %Hab2 Same vs Hab1 Diff
                    MeanFiringCellReg_30_Sh2Dh1 = cat(1,MeanFiringCellReg_30_Sh2Dh1,[mean(prev4b.catActivityHab2(pass.pass_Sh2_Dh1(:,1),:),2),mean(catActivityHab1(pass.pass_Sh2_Dh1(:,2),:),2)]);
                    MeanFiringCellReg_30_Snh2Dh1 = cat(1,MeanFiringCellReg_30_Snh2Dh1,[mean(prev4b.catActivityHab2(pass.pass_Snh2_Dh1(:,1),:),2),mean(catActivityHab1(pass.pass_Snh2_Dh1(:,2),:),2)]);
                    MeanFiringCellReg_30_Sh2Dnh1 = cat(1,MeanFiringCellReg_30_Sh2Dnh1,[mean(prev4b.catActivityHab2(pass.pass_Sh2_Dnh1(:,1),:),2),mean(catActivityHab1(pass.pass_Sh2_Dnh1(:,2),:),2)]);
                    MeanFiringCellReg_30_Snh2Dnh1 = cat(1,MeanFiringCellReg_30_Snh2Dnh1,[mean(prev4b.catActivityHab2(pass.pass_Snh2_Dnh1(:,1),:),2),mean(catActivityHab1(pass.pass_Snh2_Dnh1(:,2),:),2)]);
                    %Trial1 Same vs Trial2 Same
                    MeanFiringCellReg_30_S1S2 = cat(1,MeanFiringCellReg_30_S1D1,[mean(prev4b.catActivity1(pass.pass_S1_S2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_S1_S2(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn1S2 = cat(1,MeanFiringCellReg_30_Sn1D1,[mean(prev4b.catActivity1(pass.pass_Sn1_S2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_Sn1_S2(:,1),:),2)]);
                    MeanFiringCellReg_30_S1Sn2 = cat(1,MeanFiringCellReg_30_S1Dn1,[mean(prev4b.catActivity1(pass.pass_S1_Sn2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_S1_Sn2(:,1),:),2)]);
                    MeanFiringCellReg_30_Sn1Sn2 = cat(1,MeanFiringCellReg_30_Sn1Dn1,[mean(prev4b.catActivity1(pass.pass_Sn1_Sn2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_Sn1_Sn2(:,1),:),2)]);
                    %Trial1 Diff vs Trial2 Diff
                    MeanFiringCellReg_30_D1D2 = cat(1,MeanFiringCellReg_30_D1D2,[mean(catActivity1(pass.pass_D1_D2(:,2),:),2),mean(catActivity2(pass.pass_D1_D2(:,2),:),2)]);
                    MeanFiringCellReg_30_Dn1D2 = cat(1,MeanFiringCellReg_30_Dn1D2,[mean(catActivity1(pass.pass_Dn1_D2(:,2),:),2),mean(catActivity2(pass.pass_Dn1_D2(:,2),:),2)]);
                    MeanFiringCellReg_30_D1Dn2 = cat(1,MeanFiringCellReg_30_D1Dn2,[mean(catActivity1(pass.pass_D1_Dn2(:,2),:),2),mean(catActivity2(pass.pass_D1_Dn2(:,2),:),2)]);
                    MeanFiringCellReg_30_Dn1Dn2 = cat(1,MeanFiringCellReg_30_Dn1Dn2,[mean(catActivity1(pass.pass_Dn1_Dn2(:,2),:),2),mean(catActivity2(pass.pass_Dn1_Dn2(:,2),:),2)]);
                    
                    if isempty(Chrono30Hab1Same)
                        Chrono30Hab1Same = chronoall{previousSesh(1),1}(regind(:,1),:);
                        Chrono301Same = chronoall{previousSesh(1),2}(regind(:,1),:);
                        Chrono30Hab2Same = chronoall{previousSesh(1),3}(regind(:,1),:);
                        Chrono302Same = chronoall{previousSesh(1),4}(regind(:,1),:);
                        
                        Chrono30Hab1Diff = chronoall{previousSesh(2),1}(regind(:,2),:);
                        Chrono301Diff = chronoall{previousSesh(2),2}(regind(:,2),:);
                        Chrono30Hab2Diff = chronoall{previousSesh(2),3}(regind(:,2),:);
                        Chrono302Diff = chronoall{previousSesh(2),4}(regind(:,2),:);
                    else
                        Chrono30Hab1Same = cat(1,Chrono30Hab1Same,chronoall{previousSesh(1),1}(regind(:,1),:));
                        Chrono301Same = cat(1,Chrono301Same,chronoall{previousSesh(1),2}(regind(:,1),:));
                        Chrono30Hab2Same = cat(1,Chrono30Hab2Same,chronoall{previousSesh(1),3}(regind(:,1),:));
                        Chrono302Same = cat(1,Chrono302Same,chronoall{previousSesh(1),4}(regind(:,1),:));
                        
                        Chrono30Hab1Diff = cat(1,Chrono30Hab1Diff,chronoall{previousSesh(2),1}(regind(:,2),:));
                        Chrono301Diff = cat(1,Chrono301Diff,chronoall{previousSesh(2),2}(regind(:,2),:));
                        Chrono30Hab2Diff = cat(1,Chrono30Hab2Diff,chronoall{previousSesh(2),3}(regind(:,2),:));
                        Chrono302Diff = cat(1,Chrono302Diff,chronoall{previousSesh(2),4}(regind(:,2),:));
                    end
                    % trial 1 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,1})
                        if i > 1 && enterzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i-1) >= secframes2 && enterzonesall{SITseshnum,1}(i) + secframes2 < length(ztrial1(1,:)) && exitzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i) >= secframes2
                            bout1_2sec_30Diff_CellReg = cat(1,bout1_2sec_30Diff_CellReg,ztrial1(gcellind(regind(:,2)),enterzonesall{SITseshnum,1}(i)-secframes2:enterzonesall{SITseshnum,1}(i)+secframes2));
                        elseif i == 1 && enterzonesall{SITseshnum,1}(1) > secframes2 && exitzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i) >= secframes2
                            bout1_2sec_30Diff_CellReg = cat(1,bout1_2sec_30Diff_CellReg,ztrial1(gcellind(regind(:,2)),enterzonesall{SITseshnum,1}(i)-secframes2:enterzonesall{SITseshnum,1}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),1})
                        if i > 1 && enterzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i-1) >= secframes2 && enterzonesall{previousSesh(1),1}(i) + secframes2 < length(prevZCal.ztrial1(1,:)) && exitzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i) >= secframes2
                            bout1_2sec_30Same_CellReg = cat(1,bout1_2sec_30Same_CellReg,prevZCal.ztrial1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),1}(i)-secframes2:enterzonesall{previousSesh(1),1}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),1}(1) > secframes2 && exitzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i) >= secframes2
                            bout1_2sec_30Same_CellReg = cat(1,bout1_2sec_30Same_CellReg,prevZCal.ztrial1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),1}(i)-secframes2:enterzonesall{previousSesh(1),1}(i)+secframes2));
                        end
                    end
                    % trial 2 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,2})
                        if i > 1 && enterzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i-1) >= secframes2 && enterzonesall{SITseshnum,2}(i) + secframes2 < length(ztrial2(1,:)) && exitzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i) >= secframes2
                            bout2_2sec_30Diff_CellReg = cat(1,bout2_2sec_30Diff_CellReg,ztrial2(gcellind(regind(:,2)),enterzonesall{SITseshnum,2}(i)-secframes2:enterzonesall{SITseshnum,2}(i)+secframes2));
                        elseif i == 1 && enterzonesall{SITseshnum,2}(1) > secframes2 && exitzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i) >= secframes2
                            bout2_2sec_30Diff_CellReg = cat(1,bout2_2sec_30Diff_CellReg,ztrial2(gcellind(regind(:,2)),enterzonesall{SITseshnum,2}(i)-secframes2:enterzonesall{SITseshnum,2}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),2})
                        if i > 1 && enterzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i-1) >= secframes2 && enterzonesall{previousSesh(1),2}(i) + secframes2 < length(prevZCal.ztrial2(1,:)) && exitzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i) >= secframes2
                            bout2_2sec_30Same_CellReg = cat(1,bout2_2sec_30Same_CellReg,prevZCal.ztrial2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),2}(i)-secframes2:enterzonesall{previousSesh(1),2}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),2}(1) > secframes2 && exitzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i) >= secframes2
                            bout2_2sec_30Same_CellReg = cat(1,bout2_2sec_30Same_CellReg,prevZCal.ztrial2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),2}(i)-secframes2:enterzonesall{previousSesh(1),2}(i)+secframes2));
                        end
                    end
                    % Hab 1 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,3})
                        if i > 1 && enterzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i-1) >= secframes2 && enterzonesall{SITseshnum,3}(i) + secframes2 < length(zhab1(1,:)) && exitzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i) >= secframes2
                            boutHab1_2sec_30Diff_CellReg = cat(1,boutHab1_2sec_30Diff_CellReg,zhab1(gcellind(regind(:,2)),enterzonesall{SITseshnum,3}(i)-secframes2:enterzonesall{SITseshnum,3}(i)+secframes2));
                        elseif i == 1 && enterzonesall{SITseshnum,3}(1) > secframes2 && exitzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i) >= secframes2
                            boutHab1_2sec_30Diff_CellReg = cat(1,boutHab1_2sec_30Diff_CellReg,zhab1(gcellind(regind(:,2)),enterzonesall{SITseshnum,3}(i)-secframes2:enterzonesall{SITseshnum,3}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),3})
                        if i > 1 && enterzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i-1) >= secframes2 && enterzonesall{previousSesh(1),3}(i) + secframes2 < length(prevZCal.zhab1(1,:)) && exitzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i) >= secframes2
                            boutHab1_2sec_30Same_CellReg = cat(1,boutHab1_2sec_30Same_CellReg,prevZCal.zhab1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),3}(i)-secframes2:enterzonesall{previousSesh(1),3}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),3}(1) > secframes2 && exitzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i) >= secframes2
                            boutHab1_2sec_30Same_CellReg = cat(1,boutHab1_2sec_30Same_CellReg,prevZCal.zhab1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),3}(i)-secframes2:enterzonesall{previousSesh(1),3}(i)+secframes2));
                        end
                    end
                    % Hab 2 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,4})
                        if i > 1 && enterzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i-1) >= secframes2 && enterzonesall{SITseshnum,4}(i) + secframes2 < length(zhab2(1,:)) && exitzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i) >= secframes2
                            boutHab2_2sec_30Diff_CellReg = cat(1,boutHab2_2sec_30Diff_CellReg,zhab2(gcellind(regind(:,2)),enterzonesall{SITseshnum,4}(i)-secframes2:enterzonesall{SITseshnum,4}(i)+secframes2));
                        elseif i == 1 && enterzonesall{SITseshnum,4}(1) > secframes2 && exitzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i) >= secframes2
                            boutHab2_2sec_30Diff_CellReg = cat(1,boutHab2_2sec_30Diff_CellReg,zhab2(gcellind(regind(:,2)),enterzonesall{SITseshnum,4}(i)-secframes2:enterzonesall{SITseshnum,4}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),4})
                        if i > 1 && enterzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i-1) >= secframes2 && enterzonesall{previousSesh(1),4}(i) + secframes2 < length(prevZCal.zhab2(1,:)) && exitzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i) >= secframes2
                            boutHab2_2sec_30Same_CellReg = cat(1,boutHab2_2sec_30Same_CellReg,prevZCal.zhab2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),4}(i)-secframes2:enterzonesall{previousSesh(1),4}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),4}(1) > secframes2 && exitzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i) >= secframes2
                            boutHab2_2sec_30Same_CellReg = cat(1,boutHab2_2sec_30Same_CellReg,prevZCal.zhab2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),4}(i)-secframes2:enterzonesall{previousSesh(1),4}(i)+secframes2));
                        end
                    end
                end
            else
                if contains(SITnovelty_mins{SITseshnum,1},'Same')
                    pass = passCellReg(regind,ind_1reg,ind_2reg,ind_h1reg,ind_h2reg,ind_1reg_not,ind_2reg_not,ind_h1reg_not,ind_h2reg_not,ind_1regprev,ind_2regprev,ind_h1regprev,ind_h2regprev,ind_1regprev_not,ind_2regprev_not,ind_h1regprev_not,ind_h2regprev_not);
                    %Trial1 Same vs Trial1 Diff
                    MeanFiringCellReg_120_S1D1 = cat(1,MeanFiringCellReg_120_S1D1,[mean(catActivity1(pass.pass_S1_D1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_S1_D1(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn1D1 = cat(1,MeanFiringCellReg_120_Sn1D1,[mean(catActivity1(pass.pass_Sn1_D1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_Sn1_D1(:,1),:),2)]);
                    MeanFiringCellReg_120_S1Dn1 = cat(1,MeanFiringCellReg_120_S1Dn1,[mean(catActivity1(pass.pass_S1_Dn1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_S1_Dn1(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn1Dn1 = cat(1,MeanFiringCellReg_120_Sn1Dn1,[mean(catActivity1(pass.pass_Sn1_Dn1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_Sn1_Dn1(:,1),:),2)]);
                    %Trial2 Same vs Trial1 Diff
                    MeanFiringCellReg_120_S2D1 = cat(1,MeanFiringCellReg_120_S2D1,[mean(catActivity2(pass.pass_S2_D1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_S2_D1(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn2D1 = cat(1,MeanFiringCellReg_120_Sn2D1,[mean(catActivity2(pass.pass_Sn2_D1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_Sn2_D1(:,1),:),2)]);
                    MeanFiringCellReg_120_S2Dn1 = cat(1,MeanFiringCellReg_120_S2Dn1,[mean(catActivity2(pass.pass_S2_Dn1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_S2_Dn1(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn2Dn1 = cat(1,MeanFiringCellReg_120_Sn2Dn1,[mean(catActivity2(pass.pass_Sn2_Dn1(:,2),:),2),mean(prev4b.catActivity1(pass.pass_Sn2_Dn1(:,1),:),2)]);
                    %Trial2 Same vs Trial2 Diff
                    MeanFiringCellReg_120_S2D2 = cat(1,MeanFiringCellReg_120_S2D2,[mean(catActivity2(pass.pass_S2_D2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_S2_D2(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn2D2 = cat(1,MeanFiringCellReg_120_Sn2D2,[mean(catActivity2(pass.pass_Sn2_D2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_Sn2_D2(:,1),:),2)]);
                    MeanFiringCellReg_120_S2Dn2 = cat(1,MeanFiringCellReg_120_S2Dn2,[mean(catActivity2(pass.pass_S2_Dn2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_S2_Dn2(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn2Dn2 = cat(1,MeanFiringCellReg_120_Sn2Dn2,[mean(catActivity2(pass.pass_Sn2_Dn2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_Sn2_Dn2(:,1),:),2)]);
                    %Trial1 Same vs Trial2 Diff
                    MeanFiringCellReg_120_S1D2 = cat(1,MeanFiringCellReg_120_S1D2,[mean(catActivity1(pass.pass_S1_D2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_S1_D2(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn1D2 = cat(1,MeanFiringCellReg_120_Sn1D2,[mean(catActivity1(pass.pass_Sn1_D2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_Sn1_D2(:,1),:),2)]);
                    MeanFiringCellReg_120_S1Dn2 = cat(1,MeanFiringCellReg_120_S1Dn2,[mean(catActivity1(pass.pass_S1_Dn2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_S1_Dn2(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn1Dn2 = cat(1,MeanFiringCellReg_120_Sn1Dn2,[mean(catActivity1(pass.pass_Sn1_Dn2(:,2),:),2),mean(prev4b.catActivity2(pass.pass_Sn1_Dn2(:,1),:),2)]);
                    %Hab1 Same vs Hab1 Diff
                    MeanFiringCellReg_120_Sh1Dh1 = cat(1,MeanFiringCellReg_120_Sh1Dh1,[mean(catActivityHab1(pass.pass_Sh1_Dh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Sh1_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh1Dh1 = cat(1,MeanFiringCellReg_120_Snh1Dh1,[mean(catActivityHab1(pass.pass_Snh1_Dh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Snh1_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_120_Sh1Dnh1 = cat(1,MeanFiringCellReg_120_Sh1Dnh1,[mean(catActivityHab1(pass.pass_Sh1_Dnh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Sh1_Dnh1(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh1Dnh1 = cat(1,MeanFiringCellReg_120_Snh1Dnh1,[mean(catActivityHab1(pass.pass_Snh1_Dnh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Snh1_Dnh1(:,1),:),2)]);
                    %Hab1 Same vs Hab2 Diff
                    MeanFiringCellReg_120_Sh1Dh2 = cat(1,MeanFiringCellReg_120_Sh1Dh2,[mean(catActivityHab1(pass.pass_Sh1_Dh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Sh1_Dh2(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh1Dh2 = cat(1,MeanFiringCellReg_120_Snh1Dh2,[mean(catActivityHab1(pass.pass_Snh1_Dh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Snh1_Dh2(:,1),:),2)]);
                    MeanFiringCellReg_120_Sh1Dnh2 = cat(1,MeanFiringCellReg_120_Sh1Dnh2,[mean(catActivityHab1(pass.pass_Sh1_Dnh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Sh1_Dnh2(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh1Dnh2 = cat(1,MeanFiringCellReg_120_Snh1Dnh2,[mean(catActivityHab1(pass.pass_Snh1_Dnh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Snh1_Dnh2(:,1),:),2)]);
                    %Hab2 Same vs Hab2 Diff
                    MeanFiringCellReg_120_Sh2Dh2 = cat(1,MeanFiringCellReg_120_Sh2Dh2,[mean(catActivityHab2(pass.pass_Sh2_Dh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Sh2_Dh2(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh2Dh2 = cat(1,MeanFiringCellReg_120_Snh2Dh2,[mean(catActivityHab2(pass.pass_Snh2_Dh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Snh2_Dh2(:,1),:),2)]);
                    MeanFiringCellReg_120_Sh2Dnh2 = cat(1,MeanFiringCellReg_120_Sh2Dnh2,[mean(catActivityHab2(pass.pass_Sh2_Dnh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Sh2_Dnh2(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh2Dnh2 = cat(1,MeanFiringCellReg_120_Snh2Dnh2,[mean(catActivityHab2(pass.pass_Snh2_Dnh2(:,2),:),2),mean(prev4b.catActivityHab2(pass.pass_Snh2_Dnh2(:,1),:),2)]);
                    %Hab2 Same vs Hab1 Diff
                    MeanFiringCellReg_120_Sh2Dh1 = cat(1,MeanFiringCellReg_120_Sh2Dh1,[mean(catActivityHab2(pass.pass_Sh2_Dh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Sh2_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh2Dh1 = cat(1,MeanFiringCellReg_120_Snh2Dh1,[mean(catActivityHab2(pass.pass_Snh2_Dh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Snh2_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_120_Sh2Dnh1 = cat(1,MeanFiringCellReg_120_Sh2Dnh1,[mean(catActivityHab2(pass.pass_Sh2_Dnh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Sh2_Dnh1(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh2Dnh1 = cat(1,MeanFiringCellReg_120_Snh2Dnh1,[mean(catActivityHab2(pass.pass_Snh2_Dnh1(:,2),:),2),mean(prev4b.catActivityHab1(pass.pass_Snh2_Dnh1(:,1),:),2)]);
                    %Trial1 Same vs Trial2 Same
                    MeanFiringCellReg_120_S1S2 = cat(1,MeanFiringCellReg_120_S1S2,[mean(catActivity1(pass.pass_S1_S2(:,2),:),2),mean(catActivity2(pass.pass_S1_S2(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn1S2 = cat(1,MeanFiringCellReg_120_Sn1S2,[mean(catActivity1(pass.pass_Sn1_S2(:,2),:),2),mean(catActivity2(pass.pass_Sn1_S2(:,2),:),2)]);
                    MeanFiringCellReg_120_S1Sn2 = cat(1,MeanFiringCellReg_120_S1Sn2,[mean(catActivity1(pass.pass_S1_Sn2(:,2),:),2),mean(catActivity2(pass.pass_S1_Sn2(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn1Sn2 = cat(1,MeanFiringCellReg_120_Sn1Sn2,[mean(catActivity1(pass.pass_Sn1_Sn2(:,2),:),2),mean(catActivity2(pass.pass_Sn1_Sn2(:,2),:),2)]);
                    %Trial1 Diff vs Trial2 Diff
                    MeanFiringCellReg_120_D1D2 = cat(1,MeanFiringCellReg_120_D1D2,[mean(prev4b.catActivity1(pass.pass_D1_D2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_D1_D2(:,1),:),2)]);
                    MeanFiringCellReg_120_Dn1D2 = cat(1,MeanFiringCellReg_120_Dn1D2,[mean(prev4b.catActivity1(pass.pass_Dn1_D2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_Dn1_D2(:,1),:),2)]);
                    MeanFiringCellReg_120_D1Dn2 = cat(1,MeanFiringCellReg_120_D1Dn2,[mean(prev4b.catActivity1(pass.pass_D1_Dn2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_D1_Dn2(:,1),:),2)]);
                    MeanFiringCellReg_120_Dn1Dn2 = cat(1,MeanFiringCellReg_120_Dn1Dn2,[mean(prev4b.catActivity1(pass.pass_Dn1_Dn2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_Dn1_Dn2(:,1),:),2)]);
                    
                    if isempty(Chrono120Hab1Same)
                        Chrono120Hab1Same = chronoall{previousSesh(2),1}(regind(:,2),:);
                        Chrono1201Same = chronoall{previousSesh(2),2}(regind(:,2),:);
                        Chrono120Hab2Same = chronoall{previousSesh(2),3}(regind(:,2),:);
                        Chrono1202Same = chronoall{previousSesh(2),4}(regind(:,2),:);
                        
                        Chrono120Hab1Diff = chronoall{previousSesh(1),1}(regind(:,1),:);
                        Chrono1201Diff = chronoall{previousSesh(1),2}(regind(:,1),:);
                        Chrono120Hab2Diff = chronoall{previousSesh(1),3}(regind(:,1),:);
                        Chrono1202Diff = chronoall{previousSesh(1),4}(regind(:,1),:);
                    else
                        Chrono120Hab1Same = cat(1,Chrono120Hab1Same,chronoall{previousSesh(2),1}(regind(:,2),:));
                        Chrono1201Same = cat(1,Chrono1201Same,chronoall{previousSesh(2),2}(regind(:,2),:));
                        Chrono120Hab2Same = cat(1,Chrono120Hab2Same,chronoall{previousSesh(2),3}(regind(:,2),:));
                        Chrono1202Same = cat(1,Chrono1202Same,chronoall{previousSesh(2),4}(regind(:,2),:));
                        
                        Chrono120Hab1Diff = cat(1,Chrono120Hab1Diff,chronoall{previousSesh(1),1}(regind(:,1),:));
                        Chrono1201Diff = cat(1,Chrono1201Diff,chronoall{previousSesh(1),2}(regind(:,1),:));
                        Chrono120Hab2Diff = cat(1,Chrono120Hab2Diff,chronoall{previousSesh(1),3}(regind(:,1),:));
                        Chrono1202Diff = cat(1,Chrono1202Diff,chronoall{previousSesh(1),4}(regind(:,1),:));
                    end
                    % trial 1 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,1})
                        if i > 1 && enterzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i-1) >= secframes2 && enterzonesall{SITseshnum,1}(i) + secframes2 < length(ztrial1(1,:)) && exitzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i) >= secframes2
                            bout1_2sec_120Same_CellReg = cat(1,bout1_2sec_120Same_CellReg,ztrial1(gcellind(regind(:,2)),enterzonesall{SITseshnum,1}(i)-secframes2:enterzonesall{SITseshnum,1}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(2),1}(1) > secframes2 && exitzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i) >= secframes2
                            bout1_2sec_120Same_CellReg = cat(1,bout1_2sec_120Same_CellReg,ztrial1(gcellind(regind(:,2)),enterzonesall{SITseshnum,1}(i)-secframes2:enterzonesall{SITseshnum,1}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),1})
                        if i > 1 && enterzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i-1) >= secframes2 && enterzonesall{previousSesh(1),1}(i) + secframes2 < length(prevZCal.ztrial1(1,:)) && exitzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i) >= secframes2
                            bout1_2sec_120Diff_CellReg = cat(1,bout1_2sec_120Diff_CellReg,prevZCal.ztrial1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),1}(i)-secframes2:enterzonesall{previousSesh(1),1}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),1}(1) > secframes2 && exitzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i) >= secframes2
                            bout1_2sec_120Diff_CellReg = cat(1,bout1_2sec_120Diff_CellReg,prevZCal.ztrial1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),1}(i)-secframes2:enterzonesall{previousSesh(1),1}(i)+secframes2));
                        end
                    end
                    % trial 2 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,2})
                        if i > 1 && enterzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i-1) >= secframes2 && enterzonesall{SITseshnum,2}(i) + secframes2 < length(ztrial2(1,:)) && exitzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i) >= secframes2
                            bout2_2sec_120Same_CellReg = cat(1,bout2_2sec_120Same_CellReg,ztrial2(gcellind(regind(:,2)),enterzonesall{SITseshnum,2}(i)-secframes2:enterzonesall{SITseshnum,2}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(2),2}(1) > secframes2 && exitzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i) >= secframes2
                            bout2_2sec_120Same_CellReg = cat(1,bout2_2sec_120Same_CellReg,ztrial2(gcellind(regind(:,2)),enterzonesall{SITseshnum,2}(i)-secframes2:enterzonesall{SITseshnum,2}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),2})
                        if i > 1 && enterzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i-1) >= secframes2 && enterzonesall{previousSesh(1),2}(i)+secframes2 < length(prevZCal.ztrial2) && exitzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i) >= secframes2
                            bout2_2sec_120Diff_CellReg = cat(1,bout2_2sec_120Diff_CellReg,prevZCal.ztrial2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),2}(i)-secframes2:enterzonesall{previousSesh(1),2}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),2}(1) > secframes2 && exitzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i) >= secframes2
                            bout2_2sec_120Diff_CellReg = cat(1,bout2_2sec_120Diff_CellReg,prevZCal.ztrial2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),2}(i)-secframes2:enterzonesall{previousSesh(1),2}(i)+secframes2));
                        end
                    end
                    % Hab 1 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,3})
                        if i > 1 && enterzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i-1) >= secframes2 && enterzonesall{SITseshnum,3}(i) + secframes2 < length(zhab1(1,:)) && exitzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i) >= secframes2
                            boutHab1_2sec_120Same_CellReg = cat(1,boutHab1_2sec_120Same_CellReg,zhab1(gcellind(regind(:,2)),enterzonesall{SITseshnum,3}(i)-secframes2:enterzonesall{SITseshnum,3}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(2),3}(1) > secframes2 && exitzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i) >= secframes2
                            boutHab1_2sec_120Same_CellReg = cat(1,boutHab1_2sec_120Same_CellReg,zhab1(gcellind(regind(:,2)),enterzonesall{SITseshnum,3}(i)-secframes2:enterzonesall{SITseshnum,3}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),3})
                        if i > 1 && enterzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i-1) >= secframes2 && enterzonesall{previousSesh(1),3}(i) + secframes2 < length(prevZCal.zhab1(1,:)) && exitzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i) >= secframes2
                            boutHab1_2sec_120Diff_CellReg = cat(1,boutHab1_2sec_120Diff_CellReg,prevZCal.zhab1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),3}(i)-secframes2:enterzonesall{previousSesh(1),3}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),3}(1) > secframes2 && exitzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i) >= secframes2
                            boutHab1_2sec_120Diff_CellReg = cat(1,boutHab1_2sec_120Diff_CellReg,prevZCal.zhab1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),3}(i)-secframes2:enterzonesall{previousSesh(1),3}(i)+secframes2));
                        end
                    end
                    % Hab 2 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,4})
                        if i > 1 && enterzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i-1) >= secframes2 && enterzonesall{SITseshnum,4}(i) + secframes2 < length(prevZCal.zhab2(1,:)) && exitzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i) >= secframes2
                            boutHab2_2sec_120Same_CellReg = cat(1,boutHab2_2sec_120Same_CellReg,zhab2(gcellind(regind(:,2)),enterzonesall{SITseshnum,4}(i)-secframes2:enterzonesall{SITseshnum,4}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(2),4}(1) > secframes2 && exitzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i) >= secframes2
                            boutHab2_2sec_120Same_CellReg = cat(1,boutHab2_2sec_120Same_CellReg,zhab2(gcellind(regind(:,2)),enterzonesall{SITseshnum,4}(i)-secframes2:enterzonesall{SITseshnum,4}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),4})
                        if i > 1 && enterzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i-1) >= secframes2 && enterzonesall{previousSesh(1),4}(i) + secframes2 < length(zhab2(1,:)) && exitzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i) >= secframes2
                            boutHab2_2sec_120Diff_CellReg = cat(1,boutHab2_2sec_120Diff_CellReg,prevZCal.zhab2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),4}(i)-secframes2:enterzonesall{previousSesh(1),4}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),4}(1) > secframes2 && exitzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i) >= secframes2
                            boutHab2_2sec_120Diff_CellReg = cat(1,boutHab2_2sec_120Diff_CellReg,prevZCal.zhab2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),4}(i)-secframes2:enterzonesall{previousSesh(1),4}(i)+secframes2));
                        end
                    end
                else
                    pass = passCellReg(regind,ind_1regprev,ind_2regprev,ind_h1regprev,ind_h2regprev,ind_1regprev_not,ind_2regprev_not,ind_h1regprev_not,ind_h2regprev_not,ind_1reg,ind_2reg,ind_h1reg,ind_h2reg,ind_1reg_not,ind_2reg_not,ind_h1reg_not,ind_h2reg_not);                    
                    %Trial1 Same vs Trial1 Diff
                    MeanFiringCellReg_120_S1D1 = cat(1,MeanFiringCellReg_120_S1D1,[mean(prev4b.catActivity1(pass.pass_S1_D1(:,1),:),2),mean(catActivity1(pass.pass_S1_D1(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn1D1 = cat(1,MeanFiringCellReg_120_Sn1D1,[mean(prev4b.catActivity1(pass.pass_Sn1_D1(:,1),:),2),mean(catActivity1(pass.pass_Sn1_D1(:,2),:),2)]);
                    MeanFiringCellReg_120_S1Dn1 = cat(1,MeanFiringCellReg_120_S1Dn1,[mean(prev4b.catActivity1(pass.pass_S1_Dn1(:,1),:),2),mean(catActivity1(pass.pass_S1_Dn1(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn1Dn1 = cat(1,MeanFiringCellReg_120_Sn1Dn1,[mean(prev4b.catActivity1(pass.pass_Sn1_Dn1(:,1),:),2),mean(catActivity1(pass.pass_Sn1_Dn1(:,2),:),2)]);
                    %Trial2 Same vs Trial1 Diff
                    MeanFiringCellReg_120_S2D1 = cat(1,MeanFiringCellReg_120_S2D1,[mean(prev4b.catActivity2(pass.pass_S2_D1(:,1),:),2),mean(catActivity1(pass.pass_S2_D1(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn2D1 = cat(1,MeanFiringCellReg_120_Sn2D1,[mean(prev4b.catActivity2(pass.pass_Sn2_D1(:,1),:),2),mean(catActivity1(pass.pass_Sn2_D1(:,2),:),2)]);
                    MeanFiringCellReg_120_S2Dn1 = cat(1,MeanFiringCellReg_120_S2Dn1,[mean(prev4b.catActivity2(pass.pass_S2_Dn1(:,1),:),2),mean(catActivity1(pass.pass_S2_Dn1(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn2Dn1 = cat(1,MeanFiringCellReg_120_Sn2Dn1,[mean(prev4b.catActivity2(pass.pass_Sn2_Dn1(:,1),:),2),mean(catActivity1(pass.pass_Sn2_Dn1(:,2),:),2)]);
                    %Trial2 Same vs Trial2 Diff
                    MeanFiringCellReg_120_S2D2 = cat(1,MeanFiringCellReg_120_S2D2,[mean(prev4b.catActivity2(pass.pass_S2_D2(:,1),:),2),mean(catActivity2(pass.pass_S2_D2(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn2D2 = cat(1,MeanFiringCellReg_120_Sn2D2,[mean(prev4b.catActivity2(pass.pass_Sn2_D2(:,1),:),2),mean(catActivity2(pass.pass_Sn2_D2(:,2),:),2)]);
                    MeanFiringCellReg_120_S2Dn2 = cat(1,MeanFiringCellReg_120_S2Dn2,[mean(prev4b.catActivity2(pass.pass_S2_Dn2(:,1),:),2),mean(catActivity2(pass.pass_S2_Dn2(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn2Dn2 = cat(1,MeanFiringCellReg_120_Sn2Dn2,[mean(prev4b.catActivity2(pass.pass_Sn2_Dn2(:,1),:),2),mean(catActivity2(pass.pass_Sn2_Dn2(:,2),:),2)]);
                    %Trial1 Same vs Trial2 Diff
                    MeanFiringCellReg_120_S1D2 = cat(1,MeanFiringCellReg_120_S1D2,[mean(prev4b.catActivity1(pass.pass_S1_D2(:,1),:),2),mean(catActivity2(pass.pass_S1_D2(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn1D2 = cat(1,MeanFiringCellReg_120_Sn1D2,[mean(prev4b.catActivity1(pass.pass_Sn1_D2(:,1),:),2),mean(catActivity2(pass.pass_Sn1_D2(:,2),:),2)]);
                    MeanFiringCellReg_120_S1Dn2 = cat(1,MeanFiringCellReg_120_S1Dn2,[mean(prev4b.catActivity1(pass.pass_S1_Dn2(:,1),:),2),mean(catActivity2(pass.pass_S1_Dn2(:,2),:),2)]);
                    MeanFiringCellReg_120_Sn1Dn2 = cat(1,MeanFiringCellReg_120_Sn1Dn2,[mean(prev4b.catActivity1(pass.pass_Sn1_Dn2(:,1),:),2),mean(catActivity2(pass.pass_Sn1_Dn2(:,2),:),2)]);
                    %Hab1 Same vs Hab1 Diff
                    MeanFiringCellReg_120_Sh1Dh1 = cat(1,MeanFiringCellReg_120_Sh1Dh1,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Dh1(:,1),:),2),mean(catActivityHab1(pass.pass_Sh1_Dh1(:,2),:),2)]);
                    MeanFiringCellReg_120_Snh1Dh1 = cat(1,MeanFiringCellReg_120_Snh1Dh1,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Dh1(:,1),:),2),mean(catActivityHab1(pass.pass_Snh1_Dh1(:,2),:),2)]);
                    MeanFiringCellReg_120_Sh1Dnh1 = cat(1,MeanFiringCellReg_120_Sh1Dnh1,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Dnh1(:,1),:),2),mean(catActivityHab1(pass.pass_Sh1_Dnh1(:,2),:),2)]);
                    MeanFiringCellReg_120_Snh1Dnh1 = cat(1,MeanFiringCellReg_120_Snh1Dnh1,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Dnh1(:,1),:),2),mean(catActivityHab1(pass.pass_Snh1_Dnh1(:,2),:),2)]);
                    %Hab1 Same vs Hab2 Diff
                    MeanFiringCellReg_120_Sh1Dh2 = cat(1,MeanFiringCellReg_120_Sh1Dh2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Dh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sh1_Dh2(:,2),:),2)]);
                    MeanFiringCellReg_120_Snh1Dh2 = cat(1,MeanFiringCellReg_120_Snh1Dh2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Dh2(:,1),:),2),mean(catActivityHab2(pass.pass_Snh1_Dh2(:,2),:),2)]);
                    MeanFiringCellReg_120_Sh1Dnh2 = cat(1,MeanFiringCellReg_120_Sh1Dnh2,[mean(prev4b.catActivityHab1(pass.pass_Sh1_Dnh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sh1_Dnh2(:,2),:),2)]);
                    MeanFiringCellReg_120_Snh1Dnh2 = cat(1,MeanFiringCellReg_120_Snh1Dnh2,[mean(prev4b.catActivityHab1(pass.pass_Snh1_Dnh2(:,1),:),2),mean(catActivityHab2(pass.pass_Snh1_Dnh2(:,2),:),2)]);
                    %Hab2 Same vs Hab2 Diff
                    MeanFiringCellReg_120_Sh2Dh2 = cat(1,MeanFiringCellReg_120_Sh2Dh2,[mean(prev4b.catActivityHab2(pass.pass_Sh2_Dh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sh2_Dh2(:,2),:),2)]);
                    MeanFiringCellReg_120_Snh2Dh2 = cat(1,MeanFiringCellReg_120_Snh2Dh2,[mean(prev4b.catActivityHab2(pass.pass_Snh2_Dh2(:,1),:),2),mean(catActivityHab2(pass.pass_Snh2_Dh2(:,2),:),2)]);
                    MeanFiringCellReg_120_Sh2Dnh2 = cat(1,MeanFiringCellReg_120_Sh2Dnh2,[mean(prev4b.catActivityHab2(pass.pass_Sh2_Dnh2(:,1),:),2),mean(catActivityHab2(pass.pass_Sh2_Dnh2(:,2),:),2)]);
                    MeanFiringCellReg_120_Snh2Dnh2 = cat(1,MeanFiringCellReg_120_Snh2Dnh2,[mean(prev4b.catActivityHab2(pass.pass_Snh2_Dnh2(:,1),:),2),mean(catActivityHab2(pass.pass_Snh2_Dnh2(:,2),:),2)]);
                    %Hab2 Same vs Hab1 Diff
                    MeanFiringCellReg_120_Sh2Dh1 = cat(1,MeanFiringCellReg_120_Sh2Dh1,[mean(prev4b.catActivityHab2(pass.pass_Sh2_Dh1(:,1),:),2),mean(catActivityHab1(pass.pass_Sh2_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh2Dh1 = cat(1,MeanFiringCellReg_120_Snh2Dh1,[mean(prev4b.catActivityHab2(pass.pass_Snh2_Dh1(:,1),:),2),mean(catActivityHab1(pass.pass_Snh2_Dh1(:,1),:),2)]);
                    MeanFiringCellReg_120_Sh2Dnh1 = cat(1,MeanFiringCellReg_120_Sh2Dnh1,[mean(prev4b.catActivityHab2(pass.pass_Sh2_Dnh1(:,1),:),2),mean(catActivityHab1(pass.pass_Sh2_Dnh1(:,1),:),2)]);
                    MeanFiringCellReg_120_Snh2Dnh1 = cat(1,MeanFiringCellReg_120_Snh2Dnh1,[mean(prev4b.catActivityHab2(pass.pass_Snh2_Dnh1(:,1),:),2),mean(catActivityHab1(pass.pass_Snh2_Dnh1(:,1),:),2)]);
                    %Trial1 Same vs Trial2 Same
                    MeanFiringCellReg_120_S1S2 = cat(1,MeanFiringCellReg_120_S1D1,[mean(prev4b.catActivity1(pass.pass_S1_S2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_S1_S2(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn1S2 = cat(1,MeanFiringCellReg_120_Sn1D1,[mean(prev4b.catActivity1(pass.pass_Sn1_S2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_Sn1_S2(:,1),:),2)]);
                    MeanFiringCellReg_120_S1Sn2 = cat(1,MeanFiringCellReg_120_S1Dn1,[mean(prev4b.catActivity1(pass.pass_S1_Sn2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_S1_Sn2(:,1),:),2)]);
                    MeanFiringCellReg_120_Sn1Sn2 = cat(1,MeanFiringCellReg_120_Sn1Dn1,[mean(prev4b.catActivity1(pass.pass_Sn1_Sn2(:,1),:),2),mean(prev4b.catActivity2(pass.pass_Sn1_Sn2(:,1),:),2)]);
                    %Trial1 Diff vs Trial2 Diff
                    MeanFiringCellReg_120_D1D2 = cat(1,MeanFiringCellReg_120_D1D2,[mean(catActivity1(pass.pass_D1_D2(:,2),:),2),mean(catActivity2(pass.pass_D1_D2(:,2),:),2)]);
                    MeanFiringCellReg_120_Dn1D2 = cat(1,MeanFiringCellReg_120_Dn1D2,[mean(catActivity1(pass.pass_Dn1_D2(:,2),:),2),mean(catActivity2(pass.pass_Dn1_D2(:,2),:),2)]);
                    MeanFiringCellReg_120_D1Dn2 = cat(1,MeanFiringCellReg_120_D1Dn2,[mean(catActivity1(pass.pass_D1_Dn2(:,2),:),2),mean(catActivity2(pass.pass_D1_Dn2(:,2),:),2)]);
                    MeanFiringCellReg_120_Dn1Dn2 = cat(1,MeanFiringCellReg_120_Dn1Dn2,[mean(catActivity1(pass.pass_Dn1_Dn2(:,2),:),2),mean(catActivity2(pass.pass_Dn1_Dn2(:,2),:),2)]);
                    
                    if isempty(Chrono120Hab1Same)
                        Chrono120Hab1Same = chronoall{previousSesh(1),1}(regind(:,1),:);
                        Chrono1201Same = chronoall{previousSesh(1),2}(regind(:,1),:);
                        Chrono120Hab2Same = chronoall{previousSesh(1),3}(regind(:,1),:);
                        Chrono1202Same = chronoall{previousSesh(1),4}(regind(:,1),:);
                        
                        Chrono120Hab1Diff = chronoall{previousSesh(2),1}(regind(:,2),:);
                        Chrono1201Diff = chronoall{previousSesh(2),2}(regind(:,2),:);
                        Chrono120Hab2Diff = chronoall{previousSesh(2),3}(regind(:,2),:);
                        Chrono1202Diff = chronoall{previousSesh(2),4}(regind(:,2),:);
                    else
                        Chrono120Hab1Same = cat(1,Chrono120Hab1Same,chronoall{previousSesh(1),1}(regind(:,1),:));
                        Chrono1201Same = cat(1,Chrono1201Same,chronoall{previousSesh(1),2}(regind(:,1),:));
                        Chrono120Hab2Same = cat(1,Chrono120Hab2Same,chronoall{previousSesh(1),3}(regind(:,1),:));
                        Chrono1202Same = cat(1,Chrono1202Same,chronoall{previousSesh(1),4}(regind(:,1),:));
                        
                        Chrono120Hab1Diff = cat(1,Chrono120Hab1Diff,chronoall{previousSesh(2),1}(regind(:,2),:));
                        Chrono1201Diff = cat(1,Chrono1201Diff,chronoall{previousSesh(2),2}(regind(:,2),:));
                        Chrono120Hab2Diff = cat(1,Chrono120Hab2Diff,chronoall{previousSesh(2),3}(regind(:,2),:));
                        Chrono1202Diff = cat(1,Chrono1202Diff,chronoall{previousSesh(2),4}(regind(:,2),:));
                    end
                    % trial 1 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,1})
                        if i > 1 && enterzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i-1) >= secframes2 && enterzonesall{SITseshnum,1}(i) + secframes2 < length(ztrial1(1,:)) && exitzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i) >= secframes2
                            bout1_2sec_120Diff_CellReg = cat(1,bout1_2sec_120Diff_CellReg,ztrial1(gcellind(regind(:,2)),enterzonesall{SITseshnum,1}(i)-secframes2:enterzonesall{SITseshnum,1}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(2),1}(1) > secframes2 && exitzonesall{SITseshnum,1}(i) - enterzonesall{SITseshnum,1}(i) >= secframes2
                            bout1_2sec_120Diff_CellReg = cat(1,bout1_2sec_120Diff_CellReg,ztrial1(gcellind(regind(:,2)),enterzonesall{SITseshnum,1}(i)-secframes2:enterzonesall{SITseshnum,1}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),1})
                        if i > 1 && enterzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i-1) >= secframes2 && enterzonesall{previousSesh(1),1}(i) + secframes2 < length(prevZCal.ztrial1(1,:)) && exitzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i) >= secframes2
                            bout1_2sec_120Same_CellReg = cat(1,bout1_2sec_120Same_CellReg,prevZCal.ztrial1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),1}(i)-secframes2:enterzonesall{previousSesh(1),1}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),1}(1) > secframes2 && exitzonesall{previousSesh(1),1}(i) - enterzonesall{previousSesh(1),1}(i) >= secframes2
                            bout1_2sec_120Same_CellReg = cat(1,bout1_2sec_120Same_CellReg,prevZCal.ztrial1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),1}(i)-secframes2:enterzonesall{previousSesh(1),1}(i)+secframes2));
                        end
                    end
                    % trial 2 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,2})
                        if i > 1 && enterzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i-1) >= secframes2 && enterzonesall{SITseshnum,2}(i) + secframes2 < length(ztrial2(1,:)) && exitzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i) >= secframes2
                            bout2_2sec_120Diff_CellReg = cat(1,bout2_2sec_120Diff_CellReg,ztrial2(gcellind(regind(:,2)),enterzonesall{SITseshnum,2}(i)-secframes2:enterzonesall{SITseshnum,2}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(2),2}(1) > secframes2 && exitzonesall{SITseshnum,2}(i) - enterzonesall{SITseshnum,2}(i) >= secframes2
                            bout2_2sec_120Diff_CellReg = cat(1,bout2_2sec_120Diff_CellReg,ztrial2(gcellind(regind(:,2)),enterzonesall{SITseshnum,2}(i)-secframes2:enterzonesall{SITseshnum,2}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),2})
                        if i > 1 && enterzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i-1) >= secframes2 && enterzonesall{previousSesh(1),2}(i) + secframes2 < length(prevZCal.ztrial2(1,:))  && exitzonesall{previousSesh(1),2}(i) - enterzonesall{previousSesh(1),2}(i) >= secframes2
                            bout2_2sec_120Same_CellReg = cat(1,bout2_2sec_120Same_CellReg,prevZCal.ztrial2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),2}(i)-secframes2:enterzonesall{previousSesh(1),2}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),2}(1) > secframes2 && exitzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i) >= secframes2
                            bout2_2sec_120Same_CellReg = cat(1,bout2_2sec_120Same_CellReg,prevZCal.ztrial2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),2}(i)-secframes2:enterzonesall{previousSesh(1),2}(i)+secframes2));
                        end
                    end
                    % Hab 1 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,3})
                        if i > 1 && enterzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i-1) >= secframes2 && enterzonesall{SITseshnum,3}(i) + secframes2 < length(zhab1(1,:)) && exitzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i) >= secframes2
                            boutHab1_2sec_120Diff_CellReg = cat(1,boutHab1_2sec_120Diff_CellReg,zhab1(gcellind(regind(:,2)),enterzonesall{SITseshnum,3}(i)-secframes2:enterzonesall{SITseshnum,3}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(2),3}(1) > secframes2 && exitzonesall{SITseshnum,3}(i) - enterzonesall{SITseshnum,3}(i) >= secframes2
                            boutHab1_2sec_120Diff_CellReg = cat(1,boutHab1_2sec_120Diff_CellReg,zhab1(gcellind(regind(:,2)),enterzonesall{SITseshnum,3}(i)-secframes2:enterzonesall{SITseshnum,3}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),3})
                        if i > 1 && enterzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i-1) >= secframes2 && enterzonesall{previousSesh(1),3}(i) + secframes2 < length(prevZCal.zhab1(1,:)) && exitzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i) >= secframes2
                            boutHab1_2sec_120Same_CellReg = cat(1,boutHab1_2sec_120Same_CellReg,prevZCal.zhab1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),3}(i)-secframes2:enterzonesall{previousSesh(1),3}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),3}(1) > secframes2 && exitzonesall{previousSesh(1),3}(i) - enterzonesall{previousSesh(1),3}(i) >= secframes2
                            boutHab1_2sec_120Same_CellReg = cat(1,boutHab1_2sec_120Same_CellReg,prevZCal.zhab1(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),3}(i)-secframes2:enterzonesall{previousSesh(1),3}(i)+secframes2));
                        end
                    end
                    % Hab 2 across days Same & Diff 30 mins
                    for i = 1 : length(enterzonesall{SITseshnum,4})
                        if i > 1 && enterzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i-1) >= secframes2 && enterzonesall{SITseshnum,4}(i) + secframes2 < length(prevZCal.zhab2(1,:)) && exitzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i) >= secframes2
                            boutHab2_2sec_120Diff_CellReg = cat(1,boutHab2_2sec_120Diff_CellReg,zhab2(gcellind(regind(:,2)),enterzonesall{SITseshnum,4}(i)-secframes2:enterzonesall{SITseshnum,4}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(2),4}(1) > secframes2 && exitzonesall{SITseshnum,4}(i) - enterzonesall{SITseshnum,4}(i) >= secframes2
                            boutHab2_2sec_120Diff_CellReg = cat(1,boutHab2_2sec_120Diff_CellReg,zhab2(gcellind(regind(:,2)),enterzonesall{SITseshnum,4}(i)-secframes2:enterzonesall{SITseshnum,4}(i)+secframes2));
                        end
                    end
                    for i = 1 : length(enterzonesall{previousSesh(1),4})
                        if i > 1 && enterzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i-1) >= secframes2 && enterzonesall{previousSesh(1),4}(i) + secframes2 < length(prevZCal.zhab2(1,:)) && exitzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i) >= secframes2
                            boutHab2_2sec_120Same_CellReg = cat(1,boutHab2_2sec_120Same_CellReg,prevZCal.zhab2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),4}(i)-secframes2:enterzonesall{previousSesh(1),4}(i)+secframes2));
                        elseif i == 1 && enterzonesall{previousSesh(1),4}(1) > secframes2 && exitzonesall{previousSesh(1),4}(i) - enterzonesall{previousSesh(1),4}(i) >= secframes2
                            boutHab2_2sec_120Same_CellReg = cat(1,boutHab2_2sec_120Same_CellReg,prevZCal.zhab2(prevCellinds.gcellind(regind(:,1)),enterzonesall{previousSesh(1),4}(i)-secframes2:enterzonesall{previousSesh(1),4}(i)+secframes2));
                        end
                    end
                end
            end
            
            figure
            if contains(SITnovelty_mins{SITseshnum,1},'Same')
                subplot(2,3,1)
                scatter(meanFR1_2Raw,meanFR1_1Raw)
                title('Trial1 Novel vs Trial1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRraw-0.001 maxFRraw+0.001])
                xlim([minFRraw-0.001 maxFRraw+0.001])
                subplot(2,3,2)
                scatter(meanFR2_2Raw,meanFR2_1Raw)
                title('Trial2 Novel vs Trial2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRraw-0.001 maxFRraw+0.001])
                xlim([minFRraw-0.001 maxFRraw+0.001])
                subplot(2,3,3)
                scatter(meanFRHab1_2Raw,meanFRHab1_1Raw)
                title('Habituation1 Novel vs Habituation1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRraw-0.001 maxFRraw+0.001])
                xlim([minFRraw-0.001 maxFRraw+0.001])
                subplot(2,3,4)
                scatter(meanFRHab2_2Raw,meanFRHab2_1Raw)
                title('Habituation2 Novel vs Habituation2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRraw-0.001 maxFRraw+0.001])
                xlim([minFRraw-0.001 maxFRraw+0.001])
                
                subplot(2,3,5)
                scatter(meanFR1_2Raw,meanFR2_2Raw)
                title('Trial1 Familiar vs Trial2 Familiar')
                ylabel('Trial2 Familiar')
                xlabel('Trial1 Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([maxFRraw-0.1 maxFRraw+0.1])
                xlim([maxFRraw-0.1 maxFRraw+0.1])
                subplot(2,3,6)
                scatter(meanFR1_1Raw,meanFR2_1Raw)
                title('Trial1 Novel vs Trial2 Novel')
                ylabel('Trial2 Novel')
                xlabel('Trial1 Novel')
                hold on
                plot([-1,1],[-1,1])
                ylim([maxFRraw-0.1 maxFRraw+0.1])
                xlim([maxFRraw-0.1 maxFRraw+0.1])
            else
                subplot(2,3,1)
                scatter(meanFR1_1Raw,meanFR1_2Raw)
                title('Trial1 Novel vs Trial1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRraw-0.001 maxFRraw+0.001])
                xlim([minFRraw-0.001 maxFRraw+0.001])
                subplot(2,3,2)
                scatter(meanFR2_1Raw,meanFR2_2Raw)
                title('Trial2 Novel vs Trial2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRraw-0.001 maxFRraw+0.001])
                xlim([minFRraw-0.001 maxFRraw+0.001])
                subplot(2,3,3)
                scatter(meanFRHab1_1Raw,meanFRHab1_2Raw)
                title('Habituation1 Novel vs Habituation1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRraw-0.001 maxFRraw+0.001])
                xlim([minFRraw-0.001 maxFRraw+0.001])
                subplot(2,3,4)
                scatter(meanFRHab2_1Raw,meanFRHab2_2Raw)
                title('Habituation2 Novel vs Habituation2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRraw-0.001 maxFRraw+0.001])
                xlim([minFRraw-0.001 maxFRraw+0.001])                
                subplot(2,3,5)
                scatter(meanFR1_1Raw,meanFR2_1Raw)
                title('Trial1 Familiar vs Trial2 Familiar')
                ylabel('Trial2 Familiar')
                xlabel('Trial1 Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([maxFRraw-0.1 maxFRraw+0.1])
                xlim([maxFRraw-0.1 maxFRraw+0.1])
                subplot(2,3,6)
                scatter(meanFR1_2Raw,meanFR2_2Raw)
                title('Trial1 Novel vs Trial2 Novel')
                ylabel('Trial2 Novel')
                xlabel('Trial1 Novel')
                hold on
                plot([-1,1],[-1,1])
                ylim([maxFRraw-0.1 maxFRraw+0.1])
                xlim([maxFRraw-0.1 maxFRraw+0.1])
            end
            sgtitle([mousename ' Cell deconvoled Firing Rate Across Days' num2str(SITnovelty_mins{SITseshnum,2}) 'mins'])
            saveas(gcf,'CellFiringRateAcrossDays_Raw')
            saveas(gcf,'CellFiringRateAcrossDays_Raw.jpeg')
            figure
            if contains(SITnovelty_mins{SITseshnum,1},'Same')
                subplot(2,3,1)
                scatter(meanFR1_2z,meanFR1_1z)
                title('Trial1 Novel vs Trial1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])
                subplot(2,3,2)
                scatter(meanFR2_2z,meanFR2_1z)
                title('Trial2 Novel vs Trial2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])
                subplot(2,3,3)
                scatter(meanFRHab1_2z,meanFRHab1_1z)
                title('Habituation1 Novel vs Habituation1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])
                subplot(2,3,4)
                scatter(meanFRHab2_2z,meanFRHab2_1z)
                title('Habituation2 Novel vs Habituation2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])                
                subplot(2,3,5)
                scatter(meanFR1_2z,meanFR2_2z)
                title('Trial1 Familiar vs Trial2 Familiar')
                ylabel('Trial2 Familiar')
                xlabel('Trial1 Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])
                subplot(2,3,6)
                scatter(meanFR1_1z,meanFR2_1z)
                title('Trial1 Novel vs Trial2 Novel')
                ylabel('Trial2 Novel')
                xlabel('Trial1 Novel')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])
            else
                subplot(2,3,1)
                scatter(meanFR1_1z,meanFR1_2z)
                title('Trial1 Novel vs Trial1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])
                subplot(2,3,2)
                scatter(meanFR2_1z,meanFR2_2z)
                title('Trial2 Novel vs Trial2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])
                subplot(2,3,3)
                scatter(meanFRHab1_1z,meanFRHab1_2z)
                title('Habituation1 Novel vs Habituation1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])
                subplot(2,3,4)
                scatter(meanFRHab2_1z,meanFRHab2_2z)
                title('Habituation2 Novel vs Habituation2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])                
                subplot(2,3,5)
                scatter(meanFR1_1z,meanFR2_1z)
                title('Trial1 Familiar vs Trial2 Familiar')
                ylabel('Trial2 Familiar')
                xlabel('Trial1 Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])
                subplot(2,3,6)
                scatter(meanFR1_2z,meanFR2_2z)
                title('Trial1 Novel vs Trial2 Novel')
                ylabel('Trial2 Novel')
                xlabel('Trial1 Novel')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz-0.1 maxFRz+0.1])
                xlim([minFRz-0.1 maxFRz+0.1])                
            end
            sgtitle([mousename ' Cell Normalized Firing Rate Across Days' num2str(SITnovelty_mins{SITseshnum,2}) 'mins'])
            saveas(gcf,'CellFiringRateAcrossDays_Normalized')
            saveas(gcf,'CellFiringRateAcrossDays_Normalized.jpeg')
            
            figure
            if contains(SITnovelty_mins{SITseshnum,1},'Same')
                subplot(2,3,1)
                scatter(meanFR1_2z_passed_CR,meanFR1_1z_passed_CR)
                title('Trial1 Novel vs Trial1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])
                subplot(2,3,2)
                scatter(meanFR2_2z_passed_CR,meanFR2_1z_passed_CR)
                title('Trial2 Novel vs Trial2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])
                subplot(2,3,3)
                scatter(meanFRHab1_2z_passed_CR,meanFRHab1_1z_passed_CR)
                title('Habituation1 Novel vs Habituation1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])
                subplot(2,3,4)
                scatter(meanFRHab2_2z_passed_CR,meanFRHab2_1z_passed_CR)
                title('Habituation2 Novel vs Habituation2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])                
                subplot(2,3,5)
                scatter(meanFR1_2z_passed_CR,meanFR2_2z_passed_CR)
                title('Trial1 Familiar vs Trial2 Familiar')
                ylabel('Trial2 Familiar')
                xlabel('Trial1 Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])
                subplot(2,3,6)
                scatter(meanFR1_1z_passed_CR,meanFR2_1z_passed_CR)
                title('Trial1 Novel vs Trial2 Novel')
                ylabel('Trial2 Novel')
                xlabel('Trial1 Novel')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])
            else
                subplot(2,3,1)
                scatter(meanFR1_1z_passed_CR,meanFR1_2z_passed_CR)
                title('Trial1 Novel vs Trial1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])
                subplot(2,3,2)
                scatter(meanFR2_1z_passed_CR,meanFR2_2z_passed_CR)
                title('Trial2 Novel vs Trial2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])
                subplot(2,3,3)
                scatter(meanFRHab1_1z_passed_CR,meanFRHab1_2z_passed_CR)
                title('Habituation1 Novel vs Habituation1 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])
                subplot(2,3,4)
                scatter(meanFRHab2_1z_passed_CR,meanFRHab2_2z_passed_CR)
                title('Habituation2 Novel vs Habituation2 Familiar')
                ylabel('Novel')
                xlabel('Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])                
                subplot(2,3,5)
                scatter(meanFR1_1z_passed_CR,meanFR2_1z_passed_CR)
                title('Trial1 Familiar vs Trial2 Familiar')
                ylabel('Trial2 Familiar')
                xlabel('Trial1 Familiar')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])
                subplot(2,3,6)
                scatter(meanFR1_2z_passed_CR,meanFR2_2z_passed_CR)
                title('Trial1 Novel vs Trial2 Novel')
                ylabel('Trial2 Novel')
                xlabel('Trial1 Novel')
                hold on
                plot([-1,1],[-1,1])
                ylim([minFRz_passed-0.1 maxFRz_passed+0.1])
                xlim([minFRz_passed-0.1 maxFRz_passed+0.1])                
            end
            sgtitle([mousename 'Passed Cell Normalized Firing Rate Across Days' num2str(SITnovelty_mins{SITseshnum,2}) 'mins'])
            saveas(gcf,'CellFiringRateAcrossDays_Normalized_passedOnly')
            saveas(gcf,'CellFiringRateAcrossDays_Normalized_passedOnly.jpeg')
            
            if contains(SITnovelty_mins{SITseshnum,1},'Same')
                if SITnovelty_mins{SITseshnum,2} == 30

                    meanFR30Same_Hab1 = cat(1,meanFR30Same_Hab1,meanFRHab1_2z);
                    meanFR30Same_Trial1 = cat(1,meanFR30Same_Trial1,meanFR1_2z);
                    meanFR30Same_Hab2 = cat(1,meanFR30Same_Hab2,meanFRHab2_2z);
                    meanFR30Same_Trial2 = cat(1,meanFR30Same_Trial2,meanFR2_2z);
                    meanFR30Same_Hab1Raw = cat(1,meanFR30Same_Hab1Raw,meanFRHab1_2Raw);
                    meanFR30Same_Trial1Raw = cat(1,meanFR30Same_Trial1Raw,meanFR1_2Raw);
                    meanFR30Same_Hab2Raw = cat(1,meanFR30Same_Hab2Raw,meanFRHab2_2Raw);
                    meanFR30Same_Trial2Raw = cat(1,meanFR30Same_Trial2Raw,meanFR2_2Raw);
                    
                    meanFR30Diff_Hab1 = cat(1,meanFR30Diff_Hab1,meanFRHab1_1z);
                    meanFR30Diff_Trial1 = cat(1,meanFR30Diff_Trial1,meanFR1_1z);
                    meanFR30Diff_Hab2 = cat(1,meanFR30Diff_Hab2,meanFRHab2_1z);
                    meanFR30Diff_Trial2 = cat(1,meanFR30Diff_Trial2,meanFR2_1z);
                    meanFR30Diff_Hab1Raw = cat(1,meanFR30Diff_Hab1Raw,meanFRHab1_1Raw);
                    meanFR30Diff_Trial1Raw = cat(1,meanFR30Diff_Trial1Raw,meanFR1_1Raw);
                    meanFR30Diff_Hab2Raw = cat(1,meanFR30Diff_Hab2Raw,meanFRHab2_1Raw);
                    meanFR30Diff_Trial2Raw = cat(1,meanFR30Diff_Trial2Raw,meanFR2_1Raw);
                    
                    meanFR30Same_Hab1_passed = cat(1,meanFR30Same_Hab1_passed,meanFRHab1_2z_passed_CR);
                    meanFR30Same_Trial1_passed = cat(1,meanFR30Same_Trial1_passed,meanFR1_2z_passed_CR);
                    meanFR30Same_Hab2_passed = cat(1,meanFR30Same_Hab2_passed,meanFRHab2_2z_passed_CR);
                    meanFR30Same_Trial2_passed = cat(1,meanFR30Same_Trial2_passed,meanFR2_2z_passed_CR);
                    meanFR30Diff_Hab1_passed = cat(1,meanFR30Diff_Hab1_passed,meanFRHab1_1z_passed_CR);
                    meanFR30Diff_Trial1_passed = cat(1,meanFR30Diff_Trial1_passed,meanFR1_1z_passed_CR);
                    meanFR30Diff_Hab2_passed = cat(1,meanFR30Diff_Hab2_passed,meanFRHab2_1z_passed_CR);
                    meanFR30Diff_Trial2_passed = cat(1,meanFR30Diff_Trial2_passed,meanFR2_1z_passed_CR);
                    
                    meanFR30Same_Hab1All = cat(1,meanFR30Same_Hab1All,meanFRHab1_2zAll);
                    meanFR30Same_Trial1All = cat(1,meanFR30Same_Trial1All,meanFR1_2zAll);
                    meanFR30Same_Hab2All = cat(1,meanFR30Same_Hab2All,meanFRHab2_2zAll);
                    meanFR30Same_Trial2All = cat(1,meanFR30Same_Trial2All,meanFR2_2zAll);
                    meanFR30Diff_Hab1All = cat(1,meanFR30Diff_Hab1All,meanFRHab1_1zAll);
                    meanFR30Diff_Trial1All = cat(1,meanFR30Diff_Trial1All,meanFR1_1zAll);
                    meanFR30Diff_Hab2All = cat(1,meanFR30Diff_Hab2All,meanFRHab2_1zAll);
                    meanFR30Diff_Trial2All = cat(1,meanFR30Diff_Trial2All,meanFR2_1zAll);
                    
                    meanFR30Same_Hab1All_passed = cat(1,meanFR30Same_Hab1All_passed,meanFRHab1_2zAll_passed);
                    meanFR30Same_Trial1All_passed = cat(1,meanFR30Same_Trial1All_passed,meanFR1_2zAll_passed);
                    meanFR30Same_Hab2All_passed = cat(1,meanFR30Same_Hab2All_passed,meanFRHab2_2zAll_passed);
                    meanFR30Same_Trial2All_passed = cat(1,meanFR30Same_Trial2All_passed,meanFR2_2zAll_passed);
                    meanFR30Diff_Hab1All_passed = cat(1,meanFR30Diff_Hab1All_passed,meanFRHab1_1zAll_passed);
                    meanFR30Diff_Trial1All_passed = cat(1,meanFR30Diff_Trial1All_passed,meanFR1_1zAll_passed);
                    meanFR30Diff_Hab2All_passed = cat(1,meanFR30Diff_Hab2All_passed,meanFRHab2_1zAll_passed);
                    meanFR30Diff_Trial2All_passed = cat(1,meanFR30Diff_Trial2All_passed,meanFR2_1zAll_passed);
                else
                    
                    meanFR120Same_Hab1 = cat(1,meanFR120Same_Hab1,meanFRHab1_2z);
                    meanFR120Same_Trial1 = cat(1,meanFR120Same_Trial1,meanFR1_2z);
                    meanFR120Same_Hab2 = cat(1,meanFR120Same_Hab2,meanFRHab2_2z);
                    meanFR120Same_Trial2 = cat(1,meanFR120Same_Trial2,meanFR2_2z);
                    meanFR120Same_Hab1Raw = cat(1,meanFR120Same_Hab1Raw,meanFRHab1_2Raw);
                    meanFR120Same_Trial1Raw = cat(1,meanFR120Same_Trial1Raw,meanFR1_2Raw);
                    meanFR120Same_Hab2Raw = cat(1,meanFR120Same_Hab2Raw,meanFRHab2_2Raw);
                    meanFR120Same_Trial2Raw = cat(1,meanFR120Same_Trial2Raw,meanFR2_2Raw);
                    
                    meanFR120Diff_Hab1 = cat(1,meanFR120Diff_Hab1,meanFRHab1_1z);
                    meanFR120Diff_Trial1 = cat(1,meanFR120Diff_Trial1,meanFR1_1z);
                    meanFR120Diff_Hab2 = cat(1,meanFR120Diff_Hab2,meanFRHab2_1z);
                    meanFR120Diff_Trial2 = cat(1,meanFR120Diff_Trial2,meanFR2_1z);
                    meanFR120Diff_Hab1Raw = cat(1,meanFR120Diff_Hab1Raw,meanFRHab1_1Raw);
                    meanFR120Diff_Trial1Raw = cat(1,meanFR120Diff_Trial1Raw,meanFR1_1Raw);
                    meanFR120Diff_Hab2Raw = cat(1,meanFR120Diff_Hab2Raw,meanFRHab2_1Raw);
                    meanFR120Diff_Trial2Raw = cat(1,meanFR120Diff_Trial2Raw,meanFR2_1Raw);                                        
                    
                    meanFR120Same_Hab1_passed = cat(1,meanFR120Same_Hab1_passed,meanFRHab1_2z_passed_CR);
                    meanFR120Same_Trial1_passed = cat(1,meanFR120Same_Trial1_passed,meanFR1_2z_passed_CR);
                    meanFR120Same_Hab2_passed = cat(1,meanFR120Same_Hab2_passed,meanFRHab2_2z_passed_CR);
                    meanFR120Same_Trial2_passed = cat(1,meanFR120Same_Trial2_passed,meanFR2_2z_passed_CR);
                    meanFR120Diff_Hab1_passed = cat(1,meanFR120Diff_Hab1_passed,meanFRHab1_1z_passed_CR);
                    meanFR120Diff_Trial1_passed = cat(1,meanFR120Diff_Trial1_passed,meanFR1_1z_passed_CR);
                    meanFR120Diff_Hab2_passed = cat(1,meanFR120Diff_Hab2_passed,meanFRHab2_1z_passed_CR);
                    meanFR120Diff_Trial2_passed = cat(1,meanFR120Diff_Trial2_passed,meanFR2_1z_passed_CR);
                    
                    meanFR120Same_Hab1All = cat(1,meanFR120Same_Hab1All,meanFRHab1_2zAll);
                    meanFR120Same_Trial1All = cat(1,meanFR120Same_Trial1All,meanFR1_2zAll);
                    meanFR120Same_Hab2All = cat(1,meanFR120Same_Hab2All,meanFRHab2_2zAll);
                    meanFR120Same_Trial2All = cat(1,meanFR120Same_Trial2All,meanFR2_2zAll);
                    meanFR120Diff_Hab1All = cat(1,meanFR120Diff_Hab1All,meanFRHab1_1zAll);
                    meanFR120Diff_Trial1All = cat(1,meanFR120Diff_Trial1All,meanFR1_1zAll);
                    meanFR120Diff_Hab2All = cat(1,meanFR120Diff_Hab2All,meanFRHab2_1zAll);
                    meanFR120Diff_Trial2All = cat(1,meanFR120Diff_Trial2All,meanFR2_1zAll);
                    
                    meanFR120Same_Hab1All_passed = cat(1,meanFR120Same_Hab1All_passed,meanFRHab1_2zAll_passed);
                    meanFR120Same_Trial1All_passed = cat(1,meanFR120Same_Trial1All_passed,meanFR1_2zAll_passed);
                    meanFR120Same_Hab2All_passed = cat(1,meanFR120Same_Hab2All_passed,meanFRHab2_2zAll_passed);
                    meanFR120Same_Trial2All_passed = cat(1,meanFR120Same_Trial2All_passed,meanFR2_2zAll_passed);
                    meanFR120Diff_Hab1All_passed = cat(1,meanFR120Diff_Hab1All_passed,meanFRHab1_1zAll_passed);
                    meanFR120Diff_Trial1All_passed = cat(1,meanFR120Diff_Trial1All_passed,meanFR1_1zAll_passed);
                    meanFR120Diff_Hab2All_passed = cat(1,meanFR120Diff_Hab2All_passed,meanFRHab2_1zAll_passed);
                    meanFR120Diff_Trial2All_passed = cat(1,meanFR120Diff_Trial2All_passed,meanFR2_1zAll_passed);
                    
                end
            else
                if SITnovelty_mins{SITseshnum,2} == 30

                    meanFR30Same_Hab1 = cat(1,meanFR30Same_Hab1,meanFRHab1_1z);
                    meanFR30Same_Trial1 = cat(1,meanFR30Same_Trial1,meanFR1_1z);
                    meanFR30Same_Hab2 = cat(1,meanFR30Same_Hab2,meanFRHab2_1z);
                    meanFR30Same_Trial2 = cat(1,meanFR30Same_Trial2,meanFR2_1z);
                    meanFR30Same_Hab1Raw = cat(1,meanFR30Same_Hab1Raw,meanFRHab1_1Raw);
                    meanFR30Same_Trial1Raw = cat(1,meanFR30Same_Trial1Raw,meanFR1_1Raw);
                    meanFR30Same_Hab2Raw = cat(1,meanFR30Same_Hab2Raw,meanFRHab2_1Raw);
                    meanFR30Same_Trial2Raw = cat(1,meanFR30Same_Trial2Raw,meanFR2_1Raw);
                    
                    meanFR30Diff_Hab1 = cat(1,meanFR30Diff_Hab1,meanFRHab1_2z);
                    meanFR30Diff_Trial1 = cat(1,meanFR30Diff_Trial1,meanFR1_2z);
                    meanFR30Diff_Hab2 = cat(1,meanFR30Diff_Hab2,meanFRHab2_2z);
                    meanFR30Diff_Trial2 = cat(1,meanFR30Diff_Trial2,meanFR2_2z);
                    meanFR30Diff_Hab1Raw = cat(1,meanFR30Diff_Hab1Raw,meanFRHab1_2Raw);
                    meanFR30Diff_Trial1Raw = cat(1,meanFR30Diff_Trial1Raw,meanFR1_2Raw);
                    meanFR30Diff_Hab2Raw = cat(1,meanFR30Diff_Hab2Raw,meanFRHab2_2Raw);
                    meanFR30Diff_Trial2Raw = cat(1,meanFR30Diff_Trial2Raw,meanFR2_2Raw);
                    
                    meanFR30Same_Hab1_passed = cat(1,meanFR30Same_Hab1,meanFRHab1_1z_passed_CR);
                    meanFR30Same_Trial1_passed = cat(1,meanFR30Same_Trial1,meanFR1_1z_passed_CR);
                    meanFR30Same_Hab2_passed = cat(1,meanFR30Same_Hab2,meanFRHab2_1z_passed_CR);
                    meanFR30Same_Trial2_passed = cat(1,meanFR30Same_Trial2,meanFR2_1z_passed_CR);
                    meanFR30Diff_Hab1_passed = cat(1,meanFR30Diff_Hab1,meanFRHab1_2z_passed_CR);
                    meanFR30Diff_Trial1_passed = cat(1,meanFR30Diff_Trial1,meanFR1_2z_passed_CR);
                    meanFR30Diff_Hab2_passed = cat(1,meanFR30Diff_Hab2,meanFRHab2_2z_passed_CR);
                    meanFR30Diff_Trial2_passed = cat(1,meanFR30Diff_Trial2,meanFR2_2z_passed_CR);
                                        
                    meanFR30Same_Hab1All = cat(1,meanFR30Same_Hab1All,meanFRHab1_1zAll);
                    meanFR30Same_Trial1All = cat(1,meanFR30Same_Trial1All,meanFR1_1zAll);
                    meanFR30Same_Hab2All = cat(1,meanFR30Same_Hab2All,meanFRHab2_1zAll);
                    meanFR30Same_Trial2All = cat(1,meanFR30Same_Trial2All,meanFR2_1zAll);
                    meanFR30Diff_Hab1All = cat(1,meanFR30Diff_Hab1All,meanFRHab1_2zAll);
                    meanFR30Diff_Trial1All = cat(1,meanFR30Diff_Trial1All,meanFR1_2zAll);
                    meanFR30Diff_Hab2All = cat(1,meanFR30Diff_Hab2All,meanFRHab2_2zAll);
                    meanFR30Diff_Trial2All = cat(1,meanFR30Diff_Trial2All,meanFR2_2zAll);
                    
                    meanFR30Same_Hab1All_passed = cat(1,meanFR30Same_Hab1All_passed,meanFRHab1_1zAll_passed);
                    meanFR30Same_Trial1All_passed = cat(1,meanFR30Same_Trial1All_passed,meanFR1_1zAll_passed);
                    meanFR30Same_Hab2All_passed = cat(1,meanFR30Same_Hab2All_passed,meanFRHab2_1zAll_passed);
                    meanFR30Same_Trial2All_passed = cat(1,meanFR30Same_Trial2All_passed,meanFR2_1zAll_passed);
                    meanFR30Diff_Hab1All_passed = cat(1,meanFR30Diff_Hab1All_passed,meanFRHab1_2zAll_passed);
                    meanFR30Diff_Trial1All_passed = cat(1,meanFR30Diff_Trial1All_passed,meanFR1_2zAll_passed);
                    meanFR30Diff_Hab2All_passed = cat(1,meanFR30Diff_Hab2All_passed,meanFRHab2_2zAll_passed);
                    meanFR30Diff_Trial2All_passed = cat(1,meanFR30Diff_Trial2All_passed,meanFR2_2zAll_passed);
                else
                    meanFR120Same_Hab1 = cat(1,meanFR120Same_Hab1,meanFRHab1_1z);
                    meanFR120Same_Trial1 = cat(1,meanFR120Same_Trial1,meanFR1_1z);
                    meanFR120Same_Hab2 = cat(1,meanFR120Same_Hab2,meanFRHab2_1z);
                    meanFR120Same_Trial2 = cat(1,meanFR120Same_Trial2,meanFR2_1z);
                    meanFR120Same_Hab1Raw = cat(1,meanFR120Same_Hab1Raw,meanFRHab1_1Raw);
                    meanFR120Same_Trial1Raw = cat(1,meanFR120Same_Trial1Raw,meanFR1_1Raw);
                    meanFR120Same_Hab2Raw = cat(1,meanFR120Same_Hab2Raw,meanFRHab2_1Raw);
                    meanFR120Same_Trial2Raw = cat(1,meanFR120Same_Trial2Raw,meanFR2_1Raw);
                    
                    meanFR120Diff_Hab1 = cat(1,meanFR120Diff_Hab1,meanFRHab1_2z);
                    meanFR120Diff_Trial1 = cat(1,meanFR120Diff_Trial1,meanFR1_2z);
                    meanFR120Diff_Hab2 = cat(1,meanFR120Diff_Hab2,meanFRHab2_2z);
                    meanFR120Diff_Trial2 = cat(1,meanFR120Diff_Trial2,meanFR2_2z);
                    meanFR120Diff_Hab1Raw = cat(1,meanFR120Diff_Hab1Raw,meanFRHab1_2Raw);
                    meanFR120Diff_Trial1Raw = cat(1,meanFR120Diff_Trial1Raw,meanFR1_2Raw);
                    meanFR120Diff_Hab2Raw = cat(1,meanFR120Diff_Hab2Raw,meanFRHab2_2Raw);
                    meanFR120Diff_Trial2Raw = cat(1,meanFR120Diff_Trial2Raw,meanFR2_2Raw);
                    
                    meanFR120Same_Hab1_passed = cat(1,meanFR120Same_Hab1,meanFRHab1_1z_passed_CR);
                    meanFR120Same_Trial1_passed = cat(1,meanFR120Same_Trial1,meanFR1_1z_passed_CR);
                    meanFR120Same_Hab2_passed = cat(1,meanFR120Same_Hab2,meanFRHab2_1z_passed_CR);
                    meanFR120Same_Trial2_passed = cat(1,meanFR120Same_Trial2,meanFR2_1z_passed_CR);
                    meanFR120Diff_Hab1_passed = cat(1,meanFR120Diff_Hab1,meanFRHab1_2z_passed_CR);
                    meanFR120Diff_Trial1_passed = cat(1,meanFR120Diff_Trial1,meanFR1_2z_passed_CR);
                    meanFR120Diff_Hab2_passed = cat(1,meanFR120Diff_Hab2,meanFRHab2_2z_passed_CR);
                    meanFR120Diff_Trial2_passed = cat(1,meanFR120Diff_Trial2,meanFR2_2z_passed_CR);
                    
                    meanFR120Same_Hab1All = cat(1,meanFR120Same_Hab1All,meanFRHab1_1zAll);
                    meanFR120Same_Trial1All = cat(1,meanFR120Same_Trial1All,meanFR1_1zAll);
                    meanFR120Same_Hab2All = cat(1,meanFR120Same_Hab2All,meanFRHab2_1zAll);
                    meanFR120Same_Trial2All = cat(1,meanFR120Same_Trial2All,meanFR2_1zAll);
                    meanFR120Diff_Hab1All = cat(1,meanFR120Diff_Hab1All,meanFRHab1_2zAll);
                    meanFR120Diff_Trial1All = cat(1,meanFR120Diff_Trial1All,meanFR1_2zAll);
                    meanFR120Diff_Hab2All = cat(1,meanFR120Diff_Hab2All,meanFRHab2_2zAll);
                    meanFR120Diff_Trial2All = cat(1,meanFR120Diff_Trial2All,meanFR2_2zAll);
                    
                    meanFR120Same_Hab1All_passed = cat(1,meanFR120Same_Hab1All_passed,meanFRHab1_1zAll_passed);
                    meanFR120Same_Trial1All_passed = cat(1,meanFR120Same_Trial1All_passed,meanFR1_1zAll_passed);
                    meanFR120Same_Hab2All_passed = cat(1,meanFR120Same_Hab2All_passed,meanFRHab2_1zAll_passed);
                    meanFR120Same_Trial2All_passed = cat(1,meanFR120Same_Trial2All_passed,meanFR2_1zAll_passed);
                    meanFR120Diff_Hab1All_passed = cat(1,meanFR120Diff_Hab1All_passed,meanFRHab1_2zAll_passed);
                    meanFR120Diff_Trial1All_passed = cat(1,meanFR120Diff_Trial1All_passed,meanFR1_2zAll_passed);
                    meanFR120Diff_Hab2All_passed = cat(1,meanFR120Diff_Hab2All_passed,meanFRHab2_2zAll_passed);
                    meanFR120Diff_Trial2All_passed = cat(1,meanFR120Diff_Trial2All_passed,meanFR2_2zAll_passed);
                end
            end
        end
        regmap = [];
        alignment = [];
    end
    %     end
    close all
end
%}

%% Figure 4b Across Animals (no registration) 
[~, indall30Same1] = sort(mean(ChronoAllcells30Trial1Same,2),'descend');
[~, indall30Diff1] = sort(mean(ChronoAllcells30Trial1Diff,2),'descend');
[~, indall120Same1] = sort(mean(ChronoAllcells120Trial1Same,2),'descend');
[~, indall120Diff1] = sort(mean(ChronoAllcells120Trial1Diff,2),'descend');
[~, indall30Same2] = sort(mean(ChronoAllcells30Trial2Same,2),'descend');
[~, indall30Diff2] = sort(mean(ChronoAllcells30Trial2Diff,2),'descend');
[~, indall120Same2] = sort(mean(ChronoAllcells120Trial2Same,2),'descend');
[~, indall120Diff2] = sort(mean(ChronoAllcells120Trial2Diff,2),'descend');

maxchronoallcells30Same = max(max(cat(1,ChronoAllcells30Trial1Same,ChronoAllcells30Trial2Same,ChronoAllcells30Hab1Same,ChronoAllcells30Hab2Same)));
maxchronoallcells30Diff = max(max(cat(1,ChronoAllcells30Trial1Diff,ChronoAllcells30Trial2Diff,ChronoAllcells30Hab1Diff,ChronoAllcells30Hab2Diff)));
maxchronoallcells120Same = max(max(cat(1,ChronoAllcells120Trial1Same,ChronoAllcells120Trial2Same,ChronoAllcells120Hab1Same,ChronoAllcells120Hab2Same)));
maxchronoallcells120Diff = max(max(cat(1,ChronoAllcells120Trial1Diff,ChronoAllcells120Trial2Diff,ChronoAllcells120Hab1Diff,ChronoAllcells120Hab2Diff)));

% Familiar 30 mins
figure
subplot(1,4,1)
imagesc(ChronoAllcells30Hab1Same(indall30Same1,:))
caxis([0 maxchronoallcells30Same])
colorbar
ylabel('Cell Number')
title('Hab1 Familar')
subplot(1,4,2)
imagesc(ChronoAllcells30Trial1Same(indall30Same1,:))
caxis([0 maxchronoallcells30Same])
colorbar
title('Trial1 Familar')
subplot(1,4,3)
imagesc(ChronoAllcells30Hab2Same(indall30Same1,:))
caxis([0 maxchronoallcells30Same])
colorbar
title('Hab2 Familar')
subplot(1,4,4)
imagesc(ChronoAllcells30Trial2Same(indall30Same1,:))
caxis([0 maxchronoallcells30Same])
colorbar
title('Trial2 Familar')
sgtitle(['Interaction Activity Through Time Across Animals 30 mins'])
saveas(gcf,['Figure4b30mins_sortedTrial1Same_AcrossAllMice'])
saveas(gcf,['Figure4b30mins_sortedTrial1Same_AcrossAllMice.jpeg'])

% Familiar 120 mins 
figure
subplot(1,4,1)
imagesc(ChronoAllcells120Hab1Same(indall120Same1,:))
caxis([0 maxchronoallcells120Same])
colorbar
ylabel('Cell Number')
title('Hab1 Familar')
subplot(1,4,2)
imagesc(ChronoAllcells120Trial1Same(indall120Same1,:))
caxis([0 maxchronoallcells120Same])
colorbar
title('Trial1 Familar')
subplot(1,4,3)
imagesc(ChronoAllcells120Hab2Same(indall120Same1,:))
caxis([0 maxchronoallcells120Same])
colorbar
title('Hab2 Familar')
subplot(1,4,4)
imagesc(ChronoAllcells120Trial2Same(indall120Same1,:))
caxis([0 maxchronoallcells120Same])
colorbar
title('Trial2 Familar')
sgtitle(['Interaction Activity Through Time Across Animals 120 mins'])
saveas(gcf,['Figure4b120mins_sortedTrial1Same_AcrossAllMice'])
saveas(gcf,['Figure4b120mins_sortedTrial1Same_AcrossAllMice.jpeg'])

% Novel 30 mins
figure
subplot(1,4,1)
imagesc(ChronoAllcells30Hab1Diff(indall30Diff1,:))
caxis([0 maxchronoallcells30Diff])
colorbar
ylabel('Cell Number')
title('Hab1 Novel')
subplot(1,4,2)
imagesc(ChronoAllcells30Trial1Diff(indall30Diff1,:))
caxis([0 maxchronoallcells30Diff])
colorbar
title('Trial1 Novel')
subplot(1,4,3)
imagesc(ChronoAllcells30Hab2Diff(indall30Diff1,:))
caxis([0 maxchronoallcells30Diff])
colorbar
title('Hab2 Novel')
subplot(1,4,4)
imagesc(ChronoAllcells30Trial2Diff(indall30Diff1,:))
caxis([0 maxchronoallcells30Diff])
colorbar
title('Trial2 Novel')
sgtitle(['Interaction Activity Through Time Across Animals 30 mins'])
saveas(gcf,['Figure4b30mins_sortedTrial1Novel_AcrossAllMice'])
saveas(gcf,['Figure4b30mins_sortedTrial1Novel_AcrossAllMice.jpeg'])

% Novel 120 mins 
figure
subplot(1,4,1)
imagesc(ChronoAllcells120Hab1Diff(indall120Diff1,:))
caxis([0 maxchronoallcells120Diff])
colorbar
ylabel('Cell Number')
title('Hab1 Novel')
subplot(1,4,2)
imagesc(ChronoAllcells120Trial1Diff(indall120Diff1,:))
caxis([0 maxchronoallcells120Diff])
colorbar
title('Trial1 Novel')
subplot(1,4,3)
imagesc(ChronoAllcells120Hab2Diff(indall120Diff1,:))
caxis([0 maxchronoallcells120Diff])
colorbar
title('Hab2 Novel')
subplot(1,4,4)
imagesc(ChronoAllcells120Trial2Diff(indall120Diff1,:))
caxis([0 maxchronoallcells120Diff])
colorbar
title('Trial2 Novel')
sgtitle(['Interaction Activity Through Time Across Animals 120 mins'])
saveas(gcf,['Figure4b120mins_sortedTrial1Novel_AcrossAllMice'])
saveas(gcf,['Figure4b120mins_sortedTrial1Novel_AcrossAllMice.jpeg'])

% no hab figures
[~, indall30Same1] = sort(mean(ChronoAllcells30Trial1Same,2),'descend');
[~, indall30Diff1] = sort(mean(ChronoAllcells30Trial1Diff,2),'descend');
[~, indall120Same1] = sort(mean(ChronoAllcells120Trial1Same,2),'descend');
[~, indall120Diff1] = sort(mean(ChronoAllcells120Trial1Diff,2),'descend');
[~, indall30Same2] = sort(mean(ChronoAllcells30Trial2Same,2),'descend');
[~, indall30Diff2] = sort(mean(ChronoAllcells30Trial2Diff,2),'descend');
[~, indall120Same2] = sort(mean(ChronoAllcells120Trial2Same,2),'descend');
[~, indall120Diff2] = sort(mean(ChronoAllcells120Trial2Diff,2),'descend');

maxchronoallcells30Same = max(max(cat(1,ChronoAllcells30Trial1Same,ChronoAllcells30Trial2Same)));
maxchronoallcells30Diff = max(max(cat(1,ChronoAllcells30Trial1Diff,ChronoAllcells30Trial2Diff)));
maxchronoallcells120Same = max(max(cat(1,ChronoAllcells120Trial1Same,ChronoAllcells120Trial2Same)));
maxchronoallcells120Diff = max(max(cat(1,ChronoAllcells120Trial1Diff,ChronoAllcells120Trial2Diff)));

minchronoallcells30Same = min(min(cat(1,ChronoAllcells30Trial1Same,ChronoAllcells30Trial2Same)));
minchronoallcells30Diff = min(min(cat(1,ChronoAllcells30Trial1Diff,ChronoAllcells30Trial2Diff)));
minchronoallcells120Same = min(min(cat(1,ChronoAllcells120Trial1Same,ChronoAllcells120Trial2Same)));
minchronoallcells120Diff = min(min(cat(1,ChronoAllcells120Trial1Diff,ChronoAllcells120Trial2Diff)));

% Familiar 30 mins
figure
subplot(1,2,1)
imagesc(ChronoAllcells30Trial1Same(indall30Same1,:))
caxis([minchronoallcells30Same maxchronoallcells30Same])
colorbar
title('Trial1 Familar')
subplot(1,2,2)
imagesc(ChronoAllcells30Trial2Same(indall30Same1,:))
caxis([minchronoallcells30Same maxchronoallcells30Same])
colorbar
title('Trial2 Familar')
sgtitle(['Interaction Activity Through Time Across Animals 30 mins'])
saveas(gcf,['Figure4b30mins_sortedTrial1Same_AcrossAllMice_noHab'])
saveas(gcf,['Figure4b30mins_sortedTrial1Same_AcrossAllMice_noHab.jpeg'])

% Familiar 120 mins 
figure
subplot(1,2,1)
imagesc(ChronoAllcells120Trial1Same(indall120Same1,:))
caxis([minchronoallcells120Same maxchronoallcells120Same])
colorbar
title('Trial1 Familar')
subplot(1,2,2)
imagesc(ChronoAllcells120Trial2Same(indall120Same1,:))
caxis([minchronoallcells120Same maxchronoallcells120Same])
colorbar
title('Trial2 Familar')
sgtitle(['Interaction Activity Through Time Across Animals 120 mins'])
saveas(gcf,['Figure4b120mins_sortedTrial1Same_AcrossAllMice_noHab'])
saveas(gcf,['Figure4b120mins_sortedTrial1Same_AcrossAllMice_noHab.jpeg'])

% Novel 30 mins
figure
subplot(1,2,1)
imagesc(ChronoAllcells30Trial1Diff(indall30Diff1,:))
caxis([minchronoallcells30Diff maxchronoallcells30Diff])
colorbar
title('Trial1 Novel')
subplot(1,2,2)
imagesc(ChronoAllcells30Trial2Diff(indall30Diff1,:))
caxis([minchronoallcells30Diff maxchronoallcells30Diff])
colorbar
title('Trial2 Novel')
sgtitle(['Interaction Activity Through Time Across Animals 30 mins'])
saveas(gcf,['Figure4b30mins_sortedTrial1Novel_AcrossAllMice_noHab'])
saveas(gcf,['Figure4b30mins_sortedTrial1Novel_AcrossAllMice_noHab.jpeg'])

% Novel 120 mins 
figure
subplot(1,2,1)
imagesc(ChronoAllcells120Trial1Diff(indall120Diff1,:))
caxis([minchronoallcells120Diff maxchronoallcells120Diff])
colorbar
title('Trial1 Novel')
subplot(1,2,2)
imagesc(ChronoAllcells120Trial2Diff(indall120Diff1,:))
caxis([minchronoallcells120Diff maxchronoallcells120Diff])
colorbar
title('Trial2 Novel')
sgtitle(['Interaction Activity Through Time Across Animals 120 mins'])
saveas(gcf,['Figure4b120mins_sortedTrial1Novel_AcrossAllMice_noHab'])
saveas(gcf,['Figure4b120mins_sortedTrial1Novel_AcrossAllMice_noHab.jpeg'])


%% Figure 4b Across days/animals

n = 1;
m = 8;
dplots1 = 0.01;
dplots2 = 0.01;
margin1 = 0.1;
margin2 = 0.1;

[ha, pos] = tight_subplot(n,m,[dplots1 dplots2],[margin1 margin2],[margin1 margin2]);

[~, ind30Same1] = sort(mean(Chrono301Same,2),'descend');
[~, ind30Diff1] = sort(mean(Chrono301Diff,2),'descend');
[~, ind120Same1] = sort(mean(Chrono1201Same,2),'descend');
[~, ind120Diff1] = sort(mean(Chrono1201Diff,2),'descend');
[~, ind30Same2] = sort(mean(Chrono302Same,2),'descend');
[~, ind30Diff2] = sort(mean(Chrono302Diff,2),'descend');
[~, ind120Same2] = sort(mean(Chrono1202Same,2),'descend');
[~, ind120Diff2] = sort(mean(Chrono1202Diff,2),'descend');

maxchrono120 = max(max(cat(1,Chrono1201Diff,Chrono1202Diff,Chrono1201Same,Chrono1202Same,Chrono120Hab1Diff,Chrono120Hab2Diff,Chrono120Hab1Same,Chrono120Hab2Same)));
maxchrono30 = max(max(cat(1,Chrono301Diff,Chrono302Diff,Chrono301Same,Chrono302Same,Chrono30Hab1Diff,Chrono30Hab2Diff,Chrono30Hab1Same,Chrono30Hab2Same)));
c = 0;
% figure
% subplot(1,8,1)
c = c + 1;
axes(ha(c))
imagesc(Chrono120Hab1Same(ind120Same1,:))
caxis([0 maxchrono120])
colorbar
ylabel('Cell Number')
title('Hab1 Familar')
% subplot(1,8,2)
c = c + 1;
axes(ha(c))
imagesc(Chrono1201Same(ind120Same1,:))
caxis([0 maxchrono120])
colorbar
title('Trial1 Familar')
% subplot(1,8,3)
c = c + 1;
axes(ha(c))
imagesc(Chrono120Hab2Same(ind120Same1,:))
caxis([0 maxchrono120])
colorbar
title('Hab2 Familar')
% subplot(1,8,4)
c = c + 1;
axes(ha(c))
imagesc(Chrono1202Same(ind120Same1,:))
caxis([0 maxchrono120])
colorbar
title('Trial2 Familar')
% subplot(1,8,5)
c = c + 1;
axes(ha(c))
imagesc(Chrono120Hab1Diff(ind120Same1,:))
caxis([0 maxchrono120])
colorbar
ylabel('Cell Number')
title('Hab1 Novel')
% subplot(1,8,6)
c = c + 1;
axes(ha(c))
imagesc(Chrono1201Diff(ind120Same1,:))
caxis([0 maxchrono120])
colorbar
title('Trial1 Novel')
% subplot(1,8,7)
c = c + 1;
axes(ha(c))
imagesc(Chrono120Hab2Diff(ind120Same1,:))
caxis([0 maxchrono120])
colorbar
title('Hab2 Novel')
% subplot(1,8,8)
c = c + 1;
axes(ha(c))
imagesc(Chrono1202Diff(ind120Same1,:))
caxis([0 maxchrono120])
colorbar
title('Trial2 Novel')

sgtitle(['Interaction Activity Through Time Across Days and Animals 120 mins'])
saveas(gcf,['Figure4bAcrossDaysAllMice120mins_sortedTrial1Same'])
saveas(gcf,['Figure4bAcrossDaysAllMice120mins_sortedTrial1Same.jpeg'])
c = 0;
figure
% subplot(1,8,1)
c = c + 1;
axes(ha(c))
imagesc(Chrono30Hab1Same(ind30Same1,:))
caxis([0 maxchrono30])
colorbar
ylabel('Cell Number')
title('Hab1 Familar')
% subplot(1,8,2)
c = c + 1;
axes(ha(c))
imagesc(Chrono301Same(ind30Same1,:))
caxis([0 maxchrono30])
colorbar
title('Trial1 Familar')
% subplot(1,8,3)
c = c + 1;
axes(ha(c))
imagesc(Chrono30Hab2Same(ind30Same1,:))
caxis([0 maxchrono30])
colorbar
title('Hab2 Familar')
% subplot(1,8,4)
c = c + 1;
axes(ha(c))
imagesc(Chrono302Same(ind30Same1,:))
caxis([0 maxchrono30])
colorbar
title('Trial2 Familar')
% subplot(1,8,5)
c = c + 1;
axes(ha(c))
imagesc(Chrono30Hab1Diff(ind30Same1,:))
caxis([0 maxchrono30])
colorbar
ylabel('Cell Number')
title('Hab1 Novel')
% subplot(1,8,6)
c = c + 1;
axes(ha(c))
imagesc(Chrono301Diff(ind30Same1,:))
caxis([0 maxchrono30])
colorbar
title('Trial1 Novel')
% subplot(1,8,7)
c = c + 1;
axes(ha(c))
imagesc(Chrono30Hab2Diff(ind30Same1,:))
caxis([0 maxchrono30])
colorbar
title('Hab2 Novel')
% subplot(1,8,8)
c = c + 1;
axes(ha(c))
imagesc(Chrono302Diff(ind30Same1,:))
caxis([0 maxchrono30])
colorbar
title('Trial2 Novel')

sgtitle(['Interaction Activity Through Time Across Days and Animals 30 mins'])
saveas(gcf,['Figure4bAcrossDaysAllMice30mins_sortedTrial1Same'])
saveas(gcf,['Figure4bAcrossDaysAllMice30mins_sortedTrial1Same.jpeg'])

%% Figure 4c Across Animals (no registration)
maxFR30AllRaw = max(cat(1,MeanFiringAllcells30Trial1DiffRaw,MeanFiringAllcells30Trial1SameRaw,MeanFiringAllcells30Trial2DiffRaw,MeanFiringAllcells30Trial2SameRaw));
minFR30AllRaw = min(cat(1,MeanFiringAllcells30Trial1DiffRaw,MeanFiringAllcells30Trial1SameRaw,MeanFiringAllcells30Trial2DiffRaw,MeanFiringAllcells30Trial2SameRaw));
maxFR120AllRaw = max(cat(1,MeanFiringAllcells120Trial1DiffRaw,MeanFiringAllcells120Trial1SameRaw,MeanFiringAllcells120Trial2DiffRaw,MeanFiringAllcells120Trial2SameRaw));
minFR120AllRaw = min(cat(1,MeanFiringAllcells120Trial1DiffRaw,MeanFiringAllcells120Trial1SameRaw,MeanFiringAllcells120Trial2DiffRaw,MeanFiringAllcells120Trial2SameRaw));

maxFR30All = max(cat(1,MeanFiringAllcells30Trial1Diff,MeanFiringAllcells30Trial1Same,MeanFiringAllcells30Trial2Diff,MeanFiringAllcells30Trial2Same));
minFR30All = min(cat(1,MeanFiringAllcells30Trial1Diff,MeanFiringAllcells30Trial1Same,MeanFiringAllcells30Trial2Diff,MeanFiringAllcells30Trial2Same));
maxFR120All = max(cat(1,MeanFiringAllcells120Trial1Diff,MeanFiringAllcells120Trial1Same,MeanFiringAllcells120Trial2Diff,MeanFiringAllcells120Trial2Same));
minFR120All = min(cat(1,MeanFiringAllcells120Trial1Diff,MeanFiringAllcells120Trial1Same,MeanFiringAllcells120Trial2Diff,MeanFiringAllcells120Trial2Same));

maxFRall = max(cat(1,maxFR30All,maxFR120All));
minFRall = min(cat(1,minFR30All,minFR120All));

maxFRallRaw = max(cat(1,maxFR30AllRaw,maxFR120AllRaw));
minFRallRaw = min(cat(1,minFR30AllRaw,minFR120AllRaw));

figure
subplot(1,2,1)
scatter(MeanFiringAllcells30Trial2Same,MeanFiringAllcells30Trial1Same)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRall-0.1, maxFRall+0.1],[minFRall-0.1, maxFRall+0.1])
ylim([minFRall-0.1 maxFRall+0.1])
xlim([minFRall-0.1 maxFRall+0.1])
subplot(1,2,2)
scatter(MeanFiringAllcells30Trial2Diff,MeanFiringAllcells30Trial1Diff)
title('Trial2 Novel vs Trial2 Novel')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRall-0.1, maxFRall+0.1],[minFRall-0.1, maxFRall+0.1])
ylim([minFRall-0.1 maxFRall+0.1])
xlim([minFRall-0.1 maxFRall+0.1])

sgtitle('Normalized Firing Rate Across All Animals 30 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_NoReg')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_NoReg.jpeg')

figure
subplot(1,2,1)
scatter(MeanFiringAllcells120Trial2Same,MeanFiringAllcells120Trial1Same)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRall-0.1, maxFRall+0.1],[minFRall-0.1, maxFRall+0.1])
ylim([minFRall-0.1 maxFRall+0.1])
xlim([minFRall-0.1 maxFRall+0.1])
subplot(1,2,2)
scatter(MeanFiringAllcells120Trial2Diff,MeanFiringAllcells120Trial1Diff)
title('Trial2 Novel vs Trial2 Novel')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRall-0.1, maxFRall+0.1],[minFRall-0.1, maxFRall+0.1])
ylim([minFRall-0.1 maxFRall+0.1])
xlim([minFRall-0.1 maxFRall+0.1])
sgtitle('Normalized Firing Rate Across All Animals 120 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_NoReg')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_NoReg.jpeg')

figure
subplot(1,2,1)
scatter(MeanFiringAllcells30Trial2SameRaw,MeanFiringAllcells30Trial1SameRaw)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRallRaw-0.0001, maxFRallRaw+0.0001],[minFRallRaw-0.0001, maxFRallRaw+0.0001])
ylim([minFRallRaw-0.0001 maxFRallRaw+0.0001])
xlim([minFRallRaw-0.0001 maxFRallRaw+0.0001])
subplot(1,2,2)
scatter(MeanFiringAllcells30Trial2DiffRaw,MeanFiringAllcells30Trial1DiffRaw)
title('Trial2 Novel vs Trial2 Novel')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRallRaw-0.0001, maxFRallRaw+0.0001],[minFRallRaw-0.0001, maxFRallRaw+0.0001])
ylim([minFRallRaw-0.0001 maxFRallRaw+0.0001])
xlim([minFRallRaw-0.0001 maxFRallRaw+0.0001])
sgtitle('Deconvolved Firing Rate Across All Animals 30 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Deconvolved_NoReg')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Deconvolved_NoReg.jpeg')

figure
subplot(1,2,1)
scatter(MeanFiringAllcells120Trial2SameRaw,MeanFiringAllcells120Trial1SameRaw)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRallRaw-0.0001, maxFRallRaw+0.0001],[minFRallRaw-0.0001, maxFRallRaw+0.0001])
ylim([minFRallRaw-0.0001 maxFRallRaw+0.0001])
xlim([minFRallRaw-0.0001 maxFRallRaw+0.0001])
subplot(1,2,2)
scatter(MeanFiringAllcells120Trial2DiffRaw,MeanFiringAllcells120Trial1DiffRaw)
title('Trial2 Novel vs Trial2 Novel')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRallRaw-0.0001, maxFRallRaw+0.0001],[minFRallRaw-0.0001, maxFRallRaw+0.0001])
ylim([minFRallRaw-0.0001 maxFRallRaw+0.0001])
xlim([minFRallRaw-0.0001 maxFRallRaw+0.0001])
sgtitle('Deconvolved Firing Rate Across All Animals 120 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Deconvolved_NoReg')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Deconvolved_NoReg.jpeg')

%Passed Vals
%{
figure
subplot(1,2,1)
scatter(MeanFiringAllcells30Trial2Same_passed,MeanFiringAllcells30Trial1Same_passed)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRall-0.1, maxFRall+0.1],[minFRall-0.1, maxFRall+0.1])
ylim([minFRall-0.1 maxFRall+0.1])
xlim([minFRall-0.1 maxFRall+0.1])
subplot(1,2,2)
scatter(MeanFiringAllcells30Trial2Diff_passed,MeanFiringAllcells30Trial1Diff_passed)
title('Trial2 Novel vs Trial2 Novel')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRall-0.1, maxFRall+0.1],[minFRall-0.1, maxFRall+0.1])
ylim([minFRall-0.1 maxFRall+0.1])
xlim([minFRall-0.1 maxFRall+0.1])

sgtitle('Passed Normalized Firing Rate Across All Animals 30 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_NoReg_passed')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_NoReg_passed.jpeg')

figure
subplot(1,2,1)
scatter(MeanFiringAllcells120Trial2Same_passed,MeanFiringAllcells120Trial1Same_passed)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRall-0.1, maxFRall+0.1],[minFRall-0.1, maxFRall+0.1])
ylim([minFRall-0.1 maxFRall+0.1])
xlim([minFRall-0.1 maxFRall+0.1])
subplot(1,2,2)
scatter(MeanFiringAllcells120Trial2Diff_passed,MeanFiringAllcells120Trial1Diff_passed)
title('Trial2 Novel vs Trial2 Novel')
ylabel('Trial 1')
xlabel('Trial 2')
hold on
plot([minFRall-0.1, maxFRall+0.1],[minFRall-0.1, maxFRall+0.1])
ylim([minFRall-0.1 maxFRall+0.1])
xlim([minFRall-0.1 maxFRall+0.1])
sgtitle('Passed Normalized Firing Rate Across All Animals 120 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_NoReg_passed')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_NoReg_passed.jpeg')
%}

%% Figure 4c Across days/animals
maxFRraw30 = max(cat(1,meanFR30Diff_Hab1Raw,meanFR30Diff_Hab2Raw,meanFR30Diff_Trial1Raw,meanFR30Diff_Trial2Raw,meanFR30Same_Hab1Raw,meanFR30Same_Hab2Raw,meanFR30Same_Trial1Raw,meanFR30Same_Trial2Raw));
minFRraw30 = min(cat(1,meanFR30Diff_Hab1Raw,meanFR30Diff_Hab2Raw,meanFR30Diff_Trial1Raw,meanFR30Diff_Trial2Raw,meanFR30Same_Hab1Raw,meanFR30Same_Hab2Raw,meanFR30Same_Trial1Raw,meanFR30Same_Trial2Raw));
maxFRraw120 = max(cat(1,meanFR120Diff_Hab1Raw,meanFR120Diff_Hab2Raw,meanFR120Diff_Trial1Raw,meanFR120Diff_Trial2Raw,meanFR120Same_Hab1Raw,meanFR120Same_Hab2Raw,meanFR120Same_Trial1Raw,meanFR120Same_Trial2Raw));
minFRraw120 = min(cat(1,meanFR120Diff_Hab1Raw,meanFR120Diff_Hab2Raw,meanFR120Diff_Trial1Raw,meanFR120Diff_Trial2Raw,meanFR120Same_Hab1Raw,meanFR120Same_Hab2Raw,meanFR120Same_Trial1Raw,meanFR120Same_Trial2Raw));

maxFR30 = max(cat(1,meanFR30Diff_Hab1,meanFR30Diff_Hab2,meanFR30Diff_Trial1,meanFR30Diff_Trial2,meanFR30Same_Hab1,meanFR30Same_Hab2,meanFR30Same_Trial1,meanFR30Same_Trial2));
minFR30 = min(cat(1,meanFR30Diff_Hab1,meanFR30Diff_Hab2,meanFR30Diff_Trial1,meanFR30Diff_Trial2,meanFR30Same_Hab1,meanFR30Same_Hab2,meanFR30Same_Trial1,meanFR30Same_Trial2));
maxFR120 = max(cat(1,meanFR120Diff_Hab1,meanFR120Diff_Hab2,meanFR120Diff_Trial1,meanFR120Diff_Trial2,meanFR120Same_Hab1,meanFR120Same_Hab2,meanFR120Same_Trial1,meanFR120Same_Trial2));
minFR120 = min(cat(1,meanFR120Diff_Hab1,meanFR120Diff_Hab2,meanFR120Diff_Trial1,meanFR120Diff_Trial2,meanFR120Same_Hab1,meanFR120Same_Hab2,meanFR120Same_Trial1,meanFR120Same_Trial2));

maxFRrawAll = max(cat(1,maxFRraw30,maxFRraw120));
minFRrawAll = min(cat(1,minFRraw30,minFRraw120));

maxFRAll = max(cat(1,maxFR30,maxFR120));
minFRAll = min(cat(1,minFR30,minFR120));

figure
subplot(2,3,1)
scatter(meanFR30Same_Trial1,meanFR30Diff_Trial1)
title('Trial1 Novel vs Trial1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,2)
scatter(meanFR30Same_Trial2,meanFR30Diff_Trial2)
title('Trial2 Novel vs Trial2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,3)
scatter(meanFR30Same_Hab1,meanFR30Diff_Hab1)
title('Habituation1 Novel vs Habituation1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,4)
scatter(meanFR30Same_Hab2,meanFR30Diff_Hab2)
title('Habituation2 Novel vs Habituation2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,5)
scatter(meanFR30Same_Trial1,meanFR30Same_Trial2)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial2 Familiar')
xlabel('Trial1 Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,6)
scatter(meanFR30Diff_Trial1,meanFR30Diff_Trial2)
title('Trial1 Novel vs Trial2 Novel')
ylabel('Trial2 Novel')
xlabel('Trial1 Novel')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])

sgtitle('Normalized Firing Rate Across Days All Animals 30 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized.jpeg')

figure
subplot(2,3,1)
scatter(meanFR120Same_Trial1,meanFR120Diff_Trial1)
title('Trial1 Novel vs Trial1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,2)
scatter(meanFR120Same_Trial2,meanFR120Diff_Trial2)
title('Trial2 Novel vs Trial2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,3)
scatter(meanFR120Same_Hab1,meanFR120Diff_Hab1)
title('Habituation1 Novel vs Habituation1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,4)
scatter(meanFR120Same_Hab2,meanFR120Diff_Hab2)
title('Habituation2 Novel vs Habituation2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,5)
scatter(meanFR120Same_Trial1,meanFR120Same_Trial2)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial2 Familiar')
xlabel('Trial1 Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,6)
scatter(meanFR120Diff_Trial1,meanFR120Diff_Trial2)
title('Trial1 Novel vs Trial2 Novel')
ylabel('Trial2 Novel')
xlabel('Trial1 Novel')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])

sgtitle('Normalized Firing Rate Across Days All Animals 120 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized.jpeg')

figure
subplot(2,3,1)
scatter(meanFR30Same_Trial1Raw,meanFR30Diff_Trial1Raw)
title('Trial1 Novel vs Trial1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,2)
scatter(meanFR30Same_Trial2Raw,meanFR30Diff_Trial2Raw)
title('Trial2 Novel vs Trial2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,3)
scatter(meanFR30Same_Hab1Raw,meanFR30Diff_Hab1Raw)
title('Habituation1 Novel vs Habituation1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,4)
scatter(meanFR30Same_Hab2Raw,meanFR30Diff_Hab2Raw)
title('Habituation2 Novel vs Habituation2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,5)
scatter(meanFR30Same_Trial1Raw,meanFR30Same_Trial2Raw)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial2 Familiar')
xlabel('Trial1 Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,6)
scatter(meanFR30Diff_Trial1Raw,meanFR30Diff_Trial2Raw)
title('Trial1 Novel vs Trial2 Novel')
ylabel('Trial2 Novel')
xlabel('Trial1 Novel')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])

sgtitle('Deconvolved Firing Rate Across Days All Animals 30 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Deconvolved')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Deconvolved.jpeg')

figure
subplot(2,3,1)
scatter(meanFR120Same_Trial1Raw,meanFR120Diff_Trial1Raw)
title('Trial1 Novel vs Trial1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,2)
scatter(meanFR120Same_Trial2Raw,meanFR120Diff_Trial2Raw)
title('Trial2 Novel vs Trial2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,3)
scatter(meanFR120Same_Hab1Raw,meanFR120Diff_Hab1Raw)
title('Habituation1 Novel vs Habituation1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,4)
scatter(meanFR120Same_Hab2Raw,meanFR120Diff_Hab2Raw)
title('Habituation2 Novel vs Habituation2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,5)
scatter(meanFR120Same_Trial1Raw,meanFR120Same_Trial2Raw)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial2 Familiar')
xlabel('Trial1 Familiar')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
subplot(2,3,6)
scatter(meanFR120Diff_Trial1Raw,meanFR120Diff_Trial2Raw)
title('Trial1 Novel vs Trial2 Novel')
ylabel('Trial2 Novel')
xlabel('Trial1 Novel')
hold on
plot([minFRrawAll-0.0001, maxFRrawAll+0.0001],[minFRrawAll-0.0001, maxFRrawAll+0.0001])
ylim([minFRrawAll-0.0001 maxFRrawAll+0.0001])
xlim([minFRrawAll-0.0001 maxFRrawAll+0.0001])

sgtitle('Deconvolved Firing Rate Across Days All Animals 120 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Deconvolved')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Deconvolved.jpeg')

% across animals and days for cells that passed only
figure
subplot(2,3,1)
scatter(meanFR30Same_Trial1_passed,meanFR30Diff_Trial1_passed)
title('Trial1 Novel vs Trial1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,2)
scatter(meanFR30Same_Trial2_passed,meanFR30Diff_Trial2_passed)
title('Trial2 Novel vs Trial2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,3)
scatter(meanFR30Same_Hab1_passed,meanFR30Diff_Hab1_passed)
title('Habituation1 Novel vs Habituation1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,4)
scatter(meanFR30Same_Hab2_passed,meanFR30Diff_Hab2_passed)
title('Habituation2 Novel vs Habituation2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,5)
scatter(meanFR30Same_Trial1_passed,meanFR30Same_Trial2_passed)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial2 Familiar')
xlabel('Trial1 Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,6)
scatter(meanFR30Diff_Trial1_passed,meanFR30Diff_Trial2_passed)
title('Trial1 Novel vs Trial2 Novel')
ylabel('Trial2 Novel')
xlabel('Trial1 Novel')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])

sgtitle('Normalized Firing Rate For Passed Cells Across Days All Animals 30 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_passed')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_passed.jpeg')

figure
subplot(2,3,1)
scatter(meanFR120Same_Trial1_passed,meanFR120Diff_Trial1_passed)
title('Trial1 Novel vs Trial1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,2)
scatter(meanFR120Same_Trial2_passed,meanFR120Diff_Trial2_passed)
title('Trial2 Novel vs Trial2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,3)
scatter(meanFR120Same_Hab1_passed,meanFR120Diff_Hab1_passed)
title('Habituation1 Novel vs Habituation1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,4)
scatter(meanFR120Same_Hab2_passed,meanFR120Diff_Hab2_passed)
title('Habituation2 Novel vs Habituation2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,5)
scatter(meanFR120Same_Trial1_passed,meanFR120Same_Trial2_passed)
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial2 Familiar')
xlabel('Trial1 Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,6)
scatter(meanFR120Diff_Trial1_passed,meanFR120Diff_Trial2_passed)
title('Trial1 Novel vs Trial2 Novel')
ylabel('Trial2 Novel')
xlabel('Trial1 Novel')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])

sgtitle('Normalized Firing Rate For Passed Cells Across Days All Animals 120 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_passed')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_passed.jpeg')

%%%%%%%%%%%
%passed interactions NON-CellReg

maxFR30 = max(max(cat(1,MeanFiring_30_S1S2,MeanFiring_30_Sn1S2,MeanFiring_30_Snh1S2,MeanFiring_30_S1Sn2,MeanFiring_30_S1Snh2,MeanFiring_30_Sn1Snh2,...
    MeanFiring_30_Snh1Snh2, MeanFiring_30_Sn1S2,MeanFiring_30_Snh1S2,MeanFiring_30_S1Sn2,MeanFiring_30_S1Snh2,MeanFiring_30_Sn1Snh2,...
    MeanFiring_30_Sn1Sn2,MeanFiring_30_Sh1S2,MeanFiring_30_Sh1Sn2,MeanFiring_30_Snh1Sn2,MeanFiring_30_S1Sh2,...
    MeanFiring_30_S1Sh2,MeanFiring_30_Sn1Sh2,MeanFiring_30_Sh1Sh2,MeanFiring_30_Snh1Sh2,MeanFiring_30_Sh1Snh2,...
    MeanFiring_30_D1D2,MeanFiring_30_Dn1D2,MeanFiring_30_Dnh1D2,MeanFiring_30_D1Dn2,MeanFiring_30_D1Dnh2,MeanFiring_30_Dn1Dnh2,...
    MeanFiring_30_Dnh1Dnh2, MeanFiring_30_Dn1D2,MeanFiring_30_Dnh1D2,MeanFiring_30_D1Dn2,MeanFiring_30_D1Dnh2,MeanFiring_30_Dn1Dnh2,...
    MeanFiring_30_Dn1Dn2,MeanFiring_30_Dh1D2,MeanFiring_30_Dh1Dn2,MeanFiring_30_Dnh1Dn2,MeanFiring_30_D1Dh2,...
    MeanFiring_30_D1Dh2,MeanFiring_30_Dn1Dh2,MeanFiring_30_Dh1Dh2,MeanFiring_30_Dnh1Dh2,MeanFiring_30_Dh1Dnh2)));                                                                                
                
minFR30 = min(min(cat(1,MeanFiring_30_S1S2,MeanFiring_30_Sn1S2,MeanFiring_30_Snh1S2,MeanFiring_30_S1Sn2,MeanFiring_30_S1Snh2,MeanFiring_30_Sn1Snh2,...
    MeanFiring_30_Snh1Snh2, MeanFiring_30_Sn1S2,MeanFiring_30_Snh1S2,MeanFiring_30_S1Sn2,MeanFiring_30_S1Snh2,MeanFiring_30_Sn1Snh2,...
    MeanFiring_30_Sn1Sn2,MeanFiring_30_Sh1S2,MeanFiring_30_Sh1Sn2,MeanFiring_30_Snh1Sn2,MeanFiring_30_S1Sh2,...
    MeanFiring_30_S1Sh2,MeanFiring_30_Sn1Sh2,MeanFiring_30_Sh1Sh2,MeanFiring_30_Snh1Sh2,MeanFiring_30_Sh1Snh2,...
    MeanFiring_30_D1D2,MeanFiring_30_Dn1D2,MeanFiring_30_Dnh1D2,MeanFiring_30_D1Dn2,MeanFiring_30_D1Dnh2,MeanFiring_30_Dn1Dnh2,...
    MeanFiring_30_Dnh1Dnh2, MeanFiring_30_Dn1D2,MeanFiring_30_Dnh1D2,MeanFiring_30_D1Dn2,MeanFiring_30_D1Dnh2,MeanFiring_30_Dn1Dnh2,...
    MeanFiring_30_Dn1Dn2,MeanFiring_30_Dh1D2,MeanFiring_30_Dh1Dn2,MeanFiring_30_Dnh1Dn2,MeanFiring_30_D1Dh2,...
    MeanFiring_30_D1Dh2,MeanFiring_30_Dn1Dh2,MeanFiring_30_Dh1Dh2,MeanFiring_30_Dnh1Dh2,MeanFiring_30_Dh1Dnh2)));
maxFR120 = max(max(cat(1,MeanFiring_120_S1S2,MeanFiring_120_Sn1S2,MeanFiring_120_Snh1S2,MeanFiring_120_S1Sn2,MeanFiring_120_S1Snh2,MeanFiring_120_Sn1Snh2,...
    MeanFiring_120_Snh1Snh2, MeanFiring_120_Sn1S2,MeanFiring_120_Snh1S2,MeanFiring_120_S1Sn2,MeanFiring_120_S1Snh2,MeanFiring_120_Sn1Snh2,...
    MeanFiring_120_Sn1Sn2,MeanFiring_120_Sh1S2,MeanFiring_120_Sh1Sn2,MeanFiring_120_Snh1Sn2,MeanFiring_120_S1Sh2,...
    MeanFiring_120_S1Sh2,MeanFiring_120_Sn1Sh2,MeanFiring_120_Sh1Sh2,MeanFiring_120_Snh1Sh2,MeanFiring_120_Sh1Snh2,...
    MeanFiring_120_D1D2,MeanFiring_120_Dn1D2,MeanFiring_120_Dnh1D2,MeanFiring_120_D1Dn2,MeanFiring_120_D1Dnh2,MeanFiring_120_Dn1Dnh2,...
    MeanFiring_120_Dnh1Dnh2, MeanFiring_120_Dn1D2,MeanFiring_120_Dnh1D2,MeanFiring_120_D1Dn2,MeanFiring_120_D1Dnh2,MeanFiring_120_Dn1Dnh2,...
    MeanFiring_120_Dn1Dn2,MeanFiring_120_Dh1D2,MeanFiring_120_Dh1Dn2,MeanFiring_120_Dnh1Dn2,MeanFiring_120_D1Dh2,...
    MeanFiring_120_D1Dh2,MeanFiring_120_Dn1Dh2,MeanFiring_120_Dh1Dh2,MeanFiring_120_Dnh1Dh2,MeanFiring_120_Dh1Dnh2)));
minFR120 = min(min(cat(1,MeanFiring_120_S1S2,MeanFiring_120_Sn1S2,MeanFiring_120_Snh1S2,MeanFiring_120_S1Sn2,MeanFiring_120_S1Snh2,MeanFiring_120_Sn1Snh2,...
    MeanFiring_120_Snh1Snh2, MeanFiring_120_Sn1S2,MeanFiring_120_Snh1S2,MeanFiring_120_S1Sn2,MeanFiring_120_S1Snh2,MeanFiring_120_Sn1Snh2,...
    MeanFiring_120_Sn1Sn2,MeanFiring_120_Sh1S2,MeanFiring_120_Sh1Sn2,MeanFiring_120_Snh1Sn2,MeanFiring_120_S1Sh2,...
    MeanFiring_120_S1Sh2,MeanFiring_120_Sn1Sh2,MeanFiring_120_Sh1Sh2,MeanFiring_120_Snh1Sh2,MeanFiring_120_Sh1Snh2,...
    MeanFiring_120_D1D2,MeanFiring_120_Dn1D2,MeanFiring_120_Dnh1D2,MeanFiring_120_D1Dn2,MeanFiring_120_D1Dnh2,MeanFiring_120_Dn1Dnh2,...
    MeanFiring_120_Dnh1Dnh2, MeanFiring_120_Dn1D2,MeanFiring_120_Dnh1D2,MeanFiring_120_D1Dn2,MeanFiring_120_D1Dnh2,MeanFiring_120_Dn1Dnh2,...
    MeanFiring_120_Dn1Dn2,MeanFiring_120_Dh1D2,MeanFiring_120_Dh1Dn2,MeanFiring_120_Dnh1Dn2,MeanFiring_120_D1Dh2,...
    MeanFiring_120_D1Dh2,MeanFiring_120_Dn1Dh2,MeanFiring_120_Dh1Dh2,MeanFiring_120_Dnh1Dh2,MeanFiring_120_Dh1Dnh2)));

maxFR = max(cat(1,maxFR30,maxFR120));
minFR = min(cat(1,minFR30,minFR120));

%Colour coated Scatter 
% Green (g) for cells that pass for both Familiar (S) and Novel (D)
% Red (r) for cells that didn't pass for either
% Blue (b) for cells that don't pass Familiar (S) and pass Novel (D)
% Purple (m) for cells that pass Familiar (S) and don't pass Novel (D)

figure
subplot(1,2,1)
scatter(MeanFiring_30_S1S2(:,1),MeanFiring_30_S1S2(:,2),'g')
hold on
scatter(MeanFiring_30_Sn1S2(:,1),MeanFiring_30_Sn1S2(:,2),'b')
scatter(MeanFiring_30_S1Sn2(:,1),MeanFiring_30_S1Sn2(:,2),'m')
scatter(MeanFiring_30_Sn1Sn2(:,1),MeanFiring_30_Sn1Sn2(:,2),'r')
hold on
plot([minFR-0.1, maxFR+0.1],[minFR-0.1, maxFR+0.1])
ylim([minFR-0.1 maxFR+0.1])
xlim([minFR-0.1 maxFR+0.1])
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial 1')
xlabel('Trial 2')

subplot(1,2,2)
scatter(MeanFiring_30_D1D2(:,1),MeanFiring_30_D1D2(:,2),'g')
hold on
scatter(MeanFiring_30_Dn1D2(:,1),MeanFiring_30_Dn1D2(:,2),'b')
scatter(MeanFiring_30_D1Dn2(:,1),MeanFiring_30_D1Dn2(:,2),'m')
scatter(MeanFiring_30_Dn1Dn2(:,1),MeanFiring_30_Dn1Dn2(:,2),'r')
hold on
plot([minFR-0.1, maxFR+0.1],[minFR-0.1, maxFR+0.1])
ylim([minFR-0.1 maxFR+0.1])
xlim([minFR-0.1 maxFR+0.1])
title('Trial2 Novel vs Trial2 Novel')
ylabel('Trial 1')
xlabel('Trial 2')

sgtitle('Normalized Firing Rate Across All Animals 30 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_NoReg_ColourCoatedPassedInteractions')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_NoReg_ColourCoatedPassedInteractions.jpeg')

figure
subplot(1,2,1)
scatter(MeanFiring_120_S1S2(:,1),MeanFiring_120_S1S2(:,2),'g')
hold on
scatter(MeanFiring_120_Sn1S2(:,1),MeanFiring_120_Sn1S2(:,2),'b')
scatter(MeanFiring_120_S1Sn2(:,1),MeanFiring_120_S1Sn2(:,2),'m')
scatter(MeanFiring_120_Sn1Sn2(:,1),MeanFiring_120_Sn1Sn2(:,2),'r')
hold on
plot([minFR-0.1, maxFR+0.1],[minFR-0.1, maxFR+0.1])
ylim([minFR-0.1 maxFR+0.1])
xlim([minFR-0.1 maxFR+0.1])
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial 1')
xlabel('Trial 2')

subplot(1,2,2)
scatter(MeanFiring_120_D1D2(:,1),MeanFiring_120_D1D2(:,2),'g')
hold on
scatter(MeanFiring_120_Dn1D2(:,1),MeanFiring_120_Dn1D2(:,2),'b')
scatter(MeanFiring_120_D1Dn2(:,1),MeanFiring_120_D1Dn2(:,2),'m')
scatter(MeanFiring_120_Dn1Dn2(:,1),MeanFiring_120_Dn1Dn2(:,2),'r')
hold on
plot([minFR-0.1, maxFR+0.1],[minFR-0.1, maxFR+0.1])
ylim([minFR-0.1 maxFR+0.1])
xlim([minFR-0.1 maxFR+0.1])
title('Trial2 Novel vs Trial2 Novel')
ylabel('Trial 1')
xlabel('Trial 2')
sgtitle('Normalized Firing Rate Across All Animals 120 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_NoReg_ColourCoatedPassedInteractions')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_NoReg_ColourCoatedPassedInteractions.jpeg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%passed interactions CellReg

maxFR30_CR = max(max(cat(1,MeanFiringCellReg_30_S1D1,MeanFiringCellReg_30_Sn1D1,MeanFiringCellReg_30_S1Dn1,MeanFiringCellReg_30_Sn1Dn1,...
    MeanFiringCellReg_30_S2D1,MeanFiringCellReg_30_Sn2D1,MeanFiringCellReg_30_S2Dn1,MeanFiringCellReg_30_Sn2Dn1,...
    MeanFiringCellReg_30_S1D2,MeanFiringCellReg_30_Sn1D2,MeanFiringCellReg_30_S1Dn2,MeanFiringCellReg_30_Sn1Dn2,...
    MeanFiringCellReg_30_S2D2,MeanFiringCellReg_30_Sn2D2,MeanFiringCellReg_30_S2Dn2,MeanFiringCellReg_30_Sn2Dn2,...
    MeanFiringCellReg_30_Sh2Dh2,MeanFiringCellReg_30_Snh2Dh2,MeanFiringCellReg_30_Sh2Dnh2,MeanFiringCellReg_30_Snh2Dnh2,...
    MeanFiringCellReg_30_Sh1Dh1,MeanFiringCellReg_30_Snh1Dh1,MeanFiringCellReg_30_Sh1Dnh1,MeanFiringCellReg_30_Snh1Dnh1)));
minFR30_CR = min(min(cat(1,MeanFiringCellReg_30_S1D1,MeanFiringCellReg_30_Sn1D1,MeanFiringCellReg_30_S1Dn1,MeanFiringCellReg_30_Sn1Dn1,...
    MeanFiringCellReg_30_S2D1,MeanFiringCellReg_30_Sn2D1,MeanFiringCellReg_30_S2Dn1,MeanFiringCellReg_30_Sn2Dn1,...
    MeanFiringCellReg_30_S1D2,MeanFiringCellReg_30_Sn1D2,MeanFiringCellReg_30_S1Dn2,MeanFiringCellReg_30_Sn1Dn2,...
    MeanFiringCellReg_30_S2D2,MeanFiringCellReg_30_Sn2D2,MeanFiringCellReg_30_S2Dn2,MeanFiringCellReg_30_Sn2Dn2,...
    MeanFiringCellReg_30_Sh2Dh2,MeanFiringCellReg_30_Snh2Dh2,MeanFiringCellReg_30_Sh2Dnh2,MeanFiringCellReg_30_Snh2Dnh2,...
    MeanFiringCellReg_30_Sh1Dh1,MeanFiringCellReg_30_Snh1Dh1,MeanFiringCellReg_30_Sh1Dnh1,MeanFiringCellReg_30_Snh1Dnh1)));
maxFR120_CR = max(max(cat(1,MeanFiringCellReg_120_S1D1,MeanFiringCellReg_120_Sn1D1,MeanFiringCellReg_120_S1Dn1,MeanFiringCellReg_120_Sn1Dn1,...
    MeanFiringCellReg_120_S2D1,MeanFiringCellReg_120_Sn2D1,MeanFiringCellReg_120_S2Dn1,MeanFiringCellReg_120_Sn2Dn1,...
    MeanFiringCellReg_120_S1D2,MeanFiringCellReg_120_Sn1D2,MeanFiringCellReg_120_S1Dn2,MeanFiringCellReg_120_Sn1Dn2,...
    MeanFiringCellReg_120_S2D2,MeanFiringCellReg_120_Sn2D2,MeanFiringCellReg_120_S2Dn2,MeanFiringCellReg_120_Sn2Dn2,...
    MeanFiringCellReg_120_Sh2Dh2,MeanFiringCellReg_120_Snh2Dh2,MeanFiringCellReg_120_Sh2Dnh2,MeanFiringCellReg_120_Snh2Dnh2,...
    MeanFiringCellReg_120_Sh1Dh1,MeanFiringCellReg_120_Snh1Dh1,MeanFiringCellReg_120_Sh1Dnh1,MeanFiringCellReg_120_Snh1Dnh1)));
minFR120_CR = min(min(cat(1,MeanFiringCellReg_120_S1D1,MeanFiringCellReg_120_Sn1D1,MeanFiringCellReg_120_S1Dn1,MeanFiringCellReg_120_Sn1Dn1,...
    MeanFiringCellReg_120_S2D1,MeanFiringCellReg_120_Sn2D1,MeanFiringCellReg_120_S2Dn1,MeanFiringCellReg_120_Sn2Dn1,...
    MeanFiringCellReg_120_S1D2,MeanFiringCellReg_120_Sn1D2,MeanFiringCellReg_120_S1Dn2,MeanFiringCellReg_120_Sn1Dn2,...
    MeanFiringCellReg_120_S2D2,MeanFiringCellReg_120_Sn2D2,MeanFiringCellReg_120_S2Dn2,MeanFiringCellReg_120_Sn2Dn2,...
    MeanFiringCellReg_120_Sh2Dh2,MeanFiringCellReg_120_Snh2Dh2,MeanFiringCellReg_120_Sh2Dnh2,MeanFiringCellReg_120_Snh2Dnh2,...
    MeanFiringCellReg_120_Sh1Dh1,MeanFiringCellReg_120_Snh1Dh1,MeanFiringCellReg_120_Sh1Dnh1,MeanFiringCellReg_120_Snh1Dnh1)));

maxFRAll = max(cat(1,maxFR30_CR,maxFR120_CR));
minFRAll = min(cat(1,minFR30_CR,minFR120_CR));

%Colour coated Scatter 
% Green (g) for cells that pass for both Familiar (S) and Novel (D)
% Red (r) for cells that didn't pass for either
% Blue (b) for cells that don't pass Familiar (S) and pass Novel (D)
% Purple (m) for cells that pass Familiar (S) and don't pass Novel (D)

figure
subplot(2,3,1)
scatter(MeanFiringCellReg_30_S1D1(:,1),MeanFiringCellReg_30_S1D1(:,2),'g')
hold on
scatter(MeanFiringCellReg_30_Sn1D1(:,1),MeanFiringCellReg_30_Sn1D1(:,2),'b')
scatter(MeanFiringCellReg_30_S1Dn1(:,1),MeanFiringCellReg_30_S1Dn1(:,2),'m')
scatter(MeanFiringCellReg_30_Sn1Dn1(:,1),MeanFiringCellReg_30_Sn1Dn1(:,2),'r')
title('Trial1 Novel vs Trial1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,2)
scatter(MeanFiringCellReg_30_S2D2(:,1),MeanFiringCellReg_30_S2D2(:,2),'g')
hold on
scatter(MeanFiringCellReg_30_Sn2D2(:,1),MeanFiringCellReg_30_Sn2D2(:,2),'b')
scatter(MeanFiringCellReg_30_S2Dn2(:,1),MeanFiringCellReg_30_S2Dn2(:,2),'m')
scatter(MeanFiringCellReg_30_Sn2Dn2(:,1),MeanFiringCellReg_30_Sn2Dn2(:,2),'r')
title('Trial2 Novel vs Trial2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,3)
scatter(MeanFiringCellReg_30_Sh1Dh1(:,1),MeanFiringCellReg_30_Sh1Dh1(:,2),'g')
hold on
scatter(MeanFiringCellReg_30_Snh1Dh1(:,1),MeanFiringCellReg_30_Snh1Dh1(:,2),'b')
scatter(MeanFiringCellReg_30_Sh1Dnh1(:,1),MeanFiringCellReg_30_Sh1Dnh1(:,2),'m')
scatter(MeanFiringCellReg_30_Snh1Dnh1(:,1),MeanFiringCellReg_30_Snh1Dnh1(:,2),'r')
title('Habituation1 Novel vs Habituation1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,4)
scatter(MeanFiringCellReg_30_Sh2Dh2(:,1),MeanFiringCellReg_30_Sh2Dh2(:,2),'g')
hold on
scatter(MeanFiringCellReg_30_Snh2Dh2(:,1),MeanFiringCellReg_30_Snh2Dh2(:,2),'b')
scatter(MeanFiringCellReg_30_Sh2Dnh2(:,1),MeanFiringCellReg_30_Sh2Dnh2(:,2),'m')
scatter(MeanFiringCellReg_30_Snh2Dnh2(:,1),MeanFiringCellReg_30_Snh2Dnh2(:,2),'r')
title('Habituation2 Novel vs Habituation2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,5)
scatter(MeanFiringCellReg_30_S1S2(:,1),MeanFiringCellReg_30_S1S2(:,2),'g')
hold on
scatter(MeanFiringCellReg_30_Sn1S2(:,1),MeanFiringCellReg_30_Sn1S2(:,2),'b')
scatter(MeanFiringCellReg_30_S1Sn2(:,1),MeanFiringCellReg_30_S1Sn2(:,2),'m')
scatter(MeanFiringCellReg_30_Sn1Sn2(:,1),MeanFiringCellReg_30_Sn1Sn2(:,2),'r')
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial2 Familiar')
xlabel('Trial1 Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,6)
scatter(MeanFiringCellReg_30_D1D2(:,1),MeanFiringCellReg_30_D1D2(:,2),'g')
hold on
scatter(MeanFiringCellReg_30_Dn1D2(:,1),MeanFiringCellReg_30_Dn1D2(:,2),'b')
scatter(MeanFiringCellReg_30_D1Dn2(:,1),MeanFiringCellReg_30_D1Dn2(:,2),'m')
scatter(MeanFiringCellReg_30_Dn1Dn2(:,1),MeanFiringCellReg_30_Dn1Dn2(:,2),'r')
title('Trial1 Novel vs Trial2 Novel')
ylabel('Trial2 Novel')
xlabel('Trial1 Novel')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])

sgtitle('Passed Cells Normalized Firing Rate Across Days All Animals 30 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_PassedInteractions')
saveas(gcf,'CellFiringRateAcrossDaysAllMice30mins_Normalized_PassedInteractions.jpeg')
%120 mins 
figure
subplot(2,3,1)
scatter(MeanFiringCellReg_120_S1D1(:,1),MeanFiringCellReg_120_S1D1(:,2),'g')
hold on
scatter(MeanFiringCellReg_120_Sn1D1(:,1),MeanFiringCellReg_120_Sn1D1(:,2),'b')
scatter(MeanFiringCellReg_120_S1Dn1(:,1),MeanFiringCellReg_120_S1Dn1(:,2),'m')
scatter(MeanFiringCellReg_120_Sn1Dn1(:,1),MeanFiringCellReg_120_Sn1Dn1(:,2),'r')
title('Trial1 Novel vs Trial1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,2)
scatter(MeanFiringCellReg_120_S2D2(:,1),MeanFiringCellReg_120_S2D2(:,2),'g')
hold on
scatter(MeanFiringCellReg_120_Sn2D2(:,1),MeanFiringCellReg_120_Sn2D2(:,2),'b')
scatter(MeanFiringCellReg_120_S2Dn2(:,1),MeanFiringCellReg_120_S2Dn2(:,2),'m')
scatter(MeanFiringCellReg_120_Sn2Dn2(:,1),MeanFiringCellReg_120_Sn2Dn2(:,2),'r')
title('Trial2 Novel vs Trial2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,3)
scatter(MeanFiringCellReg_120_Sh1Dh1(:,1),MeanFiringCellReg_120_Sh1Dh1(:,2),'g')
hold on
scatter(MeanFiringCellReg_120_Snh1Dh1(:,1),MeanFiringCellReg_120_Snh1Dh1(:,2),'b')
scatter(MeanFiringCellReg_120_Sh1Dnh1(:,1),MeanFiringCellReg_120_Sh1Dnh1(:,2),'m')
scatter(MeanFiringCellReg_120_Snh1Dnh1(:,1),MeanFiringCellReg_120_Snh1Dnh1(:,2),'r')
title('Habituation1 Novel vs Habituation1 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,4)
scatter(MeanFiringCellReg_120_Sh2Dh2(:,1),MeanFiringCellReg_120_Sh2Dh2(:,2),'g')
hold on
scatter(MeanFiringCellReg_120_Snh2Dh2(:,1),MeanFiringCellReg_120_Snh2Dh2(:,2),'b')
scatter(MeanFiringCellReg_120_Sh2Dnh2(:,1),MeanFiringCellReg_120_Sh2Dnh2(:,2),'m')
scatter(MeanFiringCellReg_120_Snh2Dnh2(:,1),MeanFiringCellReg_120_Snh2Dnh2(:,2),'r')
title('Habituation2 Novel vs Habituation2 Familiar')
ylabel('Novel')
xlabel('Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,5)
scatter(MeanFiringCellReg_120_S1S2(:,1),MeanFiringCellReg_120_S1S2(:,2),'g')
hold on
scatter(MeanFiringCellReg_120_Sn1S2(:,1),MeanFiringCellReg_120_Sn1S2(:,2),'b')
scatter(MeanFiringCellReg_120_S1Sn2(:,1),MeanFiringCellReg_120_S1Sn2(:,2),'m')
scatter(MeanFiringCellReg_120_Sn1Sn2(:,1),MeanFiringCellReg_120_Sn1Sn2(:,2),'r')
title('Trial1 Familiar vs Trial2 Familiar')
ylabel('Trial2 Familiar')
xlabel('Trial1 Familiar')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])
subplot(2,3,6)
scatter(MeanFiringCellReg_120_D1D2(:,1),MeanFiringCellReg_120_D1D2(:,2),'g')
hold on
scatter(MeanFiringCellReg_120_Dn1D2(:,1),MeanFiringCellReg_120_Dn1D2(:,2),'b')
scatter(MeanFiringCellReg_120_D1Dn2(:,1),MeanFiringCellReg_120_D1Dn2(:,2),'m')
scatter(MeanFiringCellReg_120_Dn1Dn2(:,1),MeanFiringCellReg_120_Dn1Dn2(:,2),'r')
title('Trial1 Novel vs Trial2 Novel')
ylabel('Trial2 Novel')
xlabel('Trial1 Novel')
hold on
plot([minFRAll-0.1, maxFRAll+0.1],[minFRAll-0.1, maxFRAll+0.1])
ylim([minFRAll-0.1 maxFRAll+0.1])
xlim([minFRAll-0.1 maxFRAll+0.1])

sgtitle('Passed Cells Normalized Firing Rate Across Days All Animals 120 mins')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_PassedInteractions')
saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Normalized_PassedInteractions.jpeg')

%% Mean,std,n 
% Cell Reg
mean120Hab1_Diff = sum(meanFR120Diff_Hab1)/length(meanFR120Diff_Hab1);
std_mean120Hab1_Diff = std(meanFR120Diff_Hab1);
n_mean120Hab1_Diff = length(meanFR120Diff_Hab1);
mean120Hab1_Same = sum(meanFR120Same_Hab1)/length(meanFR120Same_Hab1);
std_mean120Hab1_Same = std(meanFR120Same_Hab1);
n_mean120Hab1_Same = length(meanFR120Same_Hab1);

mean120Hab2_Diff = sum(meanFR120Diff_Hab2)/length(meanFR120Diff_Hab2);
std_mean120Hab2_Diff = std(meanFR120Diff_Hab2);
n_mean120Hab2_Diff = length(meanFR120Diff_Hab2);
mean120Hab2_Same = sum(meanFR120Same_Hab2)/length(meanFR120Same_Hab2);
std_mean120Hab2_Same = std(meanFR120Same_Hab2);
n_mean120Hab2_Same = length(meanFR120Same_Hab2);

mean120Trial1_Diff = sum(meanFR120Diff_Trial1)/length(meanFR120Diff_Trial1);
std_mean120Trial1_Diff = std(meanFR120Diff_Trial1);
n_mean120Trial1_Diff = length(meanFR120Diff_Trial1);
mean120Trial1_Same = sum(meanFR120Same_Trial1)/length(meanFR120Same_Trial1);
std_mean120Trial1_Same = std(meanFR120Same_Trial1);
n_mean120Trial1_Same = length(meanFR120Same_Trial1);

mean120Trial2_Diff = sum(meanFR120Diff_Trial2)/length(meanFR120Diff_Trial2);
std_mean120Trial2_Diff = std(meanFR120Diff_Trial2);
n_mean120Trial2_Diff = length(meanFR120Diff_Trial2);
mean120Trial2_Same = sum(meanFR120Same_Trial2)/length(meanFR120Same_Trial2);
std_mean120Trial2_Same = std(meanFR120Same_Trial2);
n_mean120Trial2_Same = length(meanFR120Same_Trial2);

meanFRstats120mins = [mean120Hab1_Diff,std_mean120Hab1_Diff,n_mean120Hab1_Diff;...
    mean120Hab1_Same,std_mean120Hab1_Same,n_mean120Hab1_Same;...
    mean120Trial1_Diff,std_mean120Trial1_Diff,n_mean120Trial1_Diff;...
    mean120Trial1_Same,std_mean120Trial1_Same,n_mean120Trial1_Same;...
    mean120Hab2_Diff,std_mean120Hab2_Diff,n_mean120Hab2_Diff;...
    mean120Hab2_Same,std_mean120Hab2_Same,n_mean120Hab2_Same;...
    mean120Trial2_Diff,std_mean120Trial2_Diff,n_mean120Trial2_Diff;...
    mean120Trial2_Same,std_mean120Trial2_Same,n_mean120Trial2_Same;];

mean30Hab1_Diff = sum(meanFR30Diff_Hab1)/length(meanFR30Diff_Hab1);
std_mean30Hab1_Diff = std(meanFR30Diff_Hab1);
n_mean30Hab1_Diff = length(meanFR30Diff_Hab1);
mean30Hab1_Same = sum(meanFR30Same_Hab1)/length(meanFR30Same_Hab1);
std_mean30Hab1_Same = std(meanFR30Same_Hab1);
n_mean30Hab1_Same = length(meanFR30Same_Hab1);

mean30Hab2_Diff = sum(meanFR30Diff_Hab2)/length(meanFR30Diff_Hab2);
std_mean30Hab2_Diff = std(meanFR30Diff_Hab2);
n_mean30Hab2_Diff = length(meanFR30Diff_Hab2);
mean30Hab2_Same = sum(meanFR30Same_Hab2)/length(meanFR30Same_Hab2);
std_mean30Hab2_Same = std(meanFR30Same_Hab2);
n_mean30Hab2_Same = length(meanFR30Same_Hab2);

mean30Trial1_Diff = sum(meanFR30Diff_Trial1)/length(meanFR30Diff_Trial1);
std_mean30Trial1_Diff = std(meanFR30Diff_Trial1);
n_mean30Trial1_Diff = length(meanFR30Diff_Trial1);
mean30Trial1_Same = sum(meanFR30Same_Trial1)/length(meanFR30Same_Trial1);
std_mean30Trial1_Same = std(meanFR30Same_Trial1);
n_mean30Trial1_Same = length(meanFR30Same_Trial1);

mean30Trial2_Diff = sum(meanFR30Diff_Trial2)/length(meanFR30Diff_Trial2);
std_mean30Trial2_Diff = std(meanFR30Diff_Trial2);
n_mean30Trial2_Diff = length(meanFR30Diff_Trial2);
mean30Trial2_Same = sum(meanFR30Same_Trial2)/length(meanFR30Same_Trial2);
std_mean30Trial2_Same = std(meanFR30Same_Trial2);
n_mean30Trial2_Same = length(meanFR30Same_Trial2);

meanFRstats30mins = [mean30Hab1_Diff,std_mean30Hab1_Diff,n_mean30Hab1_Diff;...
    mean30Hab1_Same,std_mean30Hab1_Same,n_mean30Hab1_Same;...
    mean30Trial1_Diff,std_mean30Trial1_Diff,n_mean30Trial1_Diff;...
    mean30Trial1_Same,std_mean30Trial1_Same,n_mean30Trial1_Same;...
    mean30Hab2_Diff,std_mean30Hab2_Diff,n_mean30Hab2_Diff;...
    mean30Hab2_Same,std_mean30Hab2_Same,n_mean30Hab2_Same;...
    mean30Trial2_Diff,std_mean30Trial2_Diff,n_mean30Trial2_Diff;...
    mean30Trial2_Same,std_mean30Trial2_Same,n_mean30Trial2_Same;];

% z only cells that passed

mean120Hab1_Diff_passed = sum(meanFR120Diff_Hab1_passed)/length(meanFR120Diff_Hab1_passed);
std_mean120Hab1_Diff_passed = std(meanFR120Diff_Hab1_passed);
n_mean120Hab1_Diff_passed = length(meanFR120Diff_Hab1_passed);
mean120Hab1_Same_passed = sum(meanFR120Same_Hab1_passed)/length(meanFR120Same_Hab1_passed);
std_mean120Hab1_Same_passed = std(meanFR120Same_Hab1_passed);
n_mean120Hab1_Same_passed = length(meanFR120Same_Hab1_passed);

mean120Hab2_Diff_passed = sum(meanFR120Diff_Hab2_passed)/length(meanFR120Diff_Hab2_passed);
std_mean120Hab2_Diff_passed = std(meanFR120Diff_Hab2_passed);
n_mean120Hab2_Diff_passed = length(meanFR120Diff_Hab2_passed);
mean120Hab2_Same_passed = sum(meanFR120Same_Hab2_passed)/length(meanFR120Same_Hab2_passed);
std_mean120Hab2_Same_passed = std(meanFR120Same_Hab2_passed);
n_mean120Hab2_Same_passed = length(meanFR120Same_Hab2_passed);

mean120Trial1_Diff_passed = sum(meanFR120Diff_Trial1_passed)/length(meanFR120Diff_Trial1_passed);
std_mean120Trial1_Diff_passed = std(meanFR120Diff_Trial1_passed);
n_mean120Trial1_Diff_passed = length(meanFR120Diff_Trial1_passed);
mean120Trial1_Same_passed = sum(meanFR120Same_Trial1_passed)/length(meanFR120Same_Trial1_passed);
std_mean120Trial1_Same_passed = std(meanFR120Same_Trial1_passed);
n_mean120Trial1_Same_passed = length(meanFR120Same_Trial1_passed);

mean120Trial2_Diff_passed = sum(meanFR120Diff_Trial2_passed)/length(meanFR120Diff_Trial2_passed);
std_mean120Trial2_Diff_passed = std(meanFR120Diff_Trial2_passed);
n_mean120Trial2_Diff_passed = length(meanFR120Diff_Trial2_passed);
mean120Trial2_Same_passed = sum(meanFR120Same_Trial2_passed)/length(meanFR120Same_Trial2_passed);
std_mean120Trial2_Same_passed = std(meanFR120Same_Trial2_passed);
n_mean120Trial2_Same_passed = length(meanFR120Same_Trial2_passed);

meanFRstats120mins_passed = [mean120Hab1_Diff_passed,std_mean120Hab1_Diff_passed,n_mean120Hab1_Diff_passed;...
    mean120Hab1_Same_passed,std_mean120Hab1_Same_passed,n_mean120Hab1_Same_passed;...
    mean120Trial1_Diff_passed,std_mean120Trial1_Diff_passed,n_mean120Trial1_Diff_passed;...
    mean120Trial1_Same_passed,std_mean120Trial1_Same_passed,n_mean120Trial1_Same_passed;...
    mean120Hab2_Diff_passed,std_mean120Hab2_Diff_passed,n_mean120Hab2_Diff_passed;...
    mean120Hab2_Same_passed,std_mean120Hab2_Same_passed,n_mean120Hab2_Same_passed;...
    mean120Trial2_Diff_passed,std_mean120Trial2_Diff_passed,n_mean120Trial2_Diff_passed;...
    mean120Trial2_Same_passed,std_mean120Trial2_Same_passed,n_mean120Trial2_Same_passed;];

mean30Hab1_Diff_passed = sum(meanFR30Diff_Hab1_passed)/length(meanFR30Diff_Hab1_passed);
std_mean30Hab1_Diff_passed = std(meanFR30Diff_Hab1_passed);
n_mean30Hab1_Diff_passed = length(meanFR30Diff_Hab1_passed);
mean30Hab1_Same_passed = sum(meanFR30Same_Hab1_passed)/length(meanFR30Same_Hab1_passed);
std_mean30Hab1_Same_passed = std(meanFR30Same_Hab1_passed);
n_mean30Hab1_Same_passed = length(meanFR30Same_Hab1_passed);

mean30Hab2_Diff_passed = sum(meanFR30Diff_Hab2_passed)/length(meanFR30Diff_Hab2_passed);
std_mean30Hab2_Diff_passed = std(meanFR30Diff_Hab2_passed);
n_mean30Hab2_Diff_passed = length(meanFR30Diff_Hab2_passed);
mean30Hab2_Same_passed = sum(meanFR30Same_Hab2_passed)/length(meanFR30Same_Hab2_passed);
std_mean30Hab2_Same_passed = std(meanFR30Same_Hab2_passed);
n_mean30Hab2_Same_passed = length(meanFR30Same_Hab2_passed);

mean30Trial1_Diff_passed = sum(meanFR30Diff_Trial1_passed)/length(meanFR30Diff_Trial1_passed);
std_mean30Trial1_Diff_passed = std(meanFR30Diff_Trial1_passed);
n_mean30Trial1_Diff_passed = length(meanFR30Diff_Trial1_passed);
mean30Trial1_Same_passed = sum(meanFR30Same_Trial1_passed)/length(meanFR30Same_Trial1_passed);
std_mean30Trial1_Same_passed = std(meanFR30Same_Trial1_passed);
n_mean30Trial1_Same_passed = length(meanFR30Same_Trial1);

mean30Trial2_Diff_passed = sum(meanFR30Diff_Trial2_passed)/length(meanFR30Diff_Trial2_passed);
std_mean30Trial2_Diff_passed = std(meanFR30Diff_Trial2_passed);
n_mean30Trial2_Diff_passed = length(meanFR30Diff_Trial2_passed);
mean30Trial2_Same_passed = sum(meanFR30Same_Trial2_passed)/length(meanFR30Same_Trial2_passed);
std_mean30Trial2_Same_passed = std(meanFR30Same_Trial2_passed);
n_mean30Trial2_Same_passed = length(meanFR30Same_Trial2_passed);

meanFRstats30mins_passed = [mean30Hab1_Diff_passed,std_mean30Hab1_Diff_passed,n_mean30Hab1_Diff_passed;...
    mean30Hab1_Same_passed,std_mean30Hab1_Same_passed,n_mean30Hab1_Same_passed;...
    mean30Trial1_Diff_passed,std_mean30Trial1_Diff_passed,n_mean30Trial1_Diff_passed;...
    mean30Trial1_Same_passed,std_mean30Trial1_Same_passed,n_mean30Trial1_Same_passed;...
    mean30Hab2_Diff_passed,std_mean30Hab2_Diff_passed,n_mean30Hab2_Diff_passed;...
    mean30Hab2_Same_passed,std_mean30Hab2_Same_passed,n_mean30Hab2_Same_passed;...
    mean30Trial2_Diff_passed,std_mean30Trial2_Diff_passed,n_mean30Trial2_Diff_passed;...
    mean30Trial2_Same_passed,std_mean30Trial2_Same_passed,n_mean30Trial2_Same_passed;];

% raw
mean120Hab1_DiffRaw = sum(meanFR120Diff_Hab1Raw)/length(meanFR120Diff_Hab1Raw);
std_mean120Hab1_DiffRaw = std(meanFR120Diff_Hab1Raw);
n_mean120Hab1_DiffRaw = length(meanFR120Diff_Hab1Raw);
mean120Hab1_SameRaw = sum(meanFR120Same_Hab1Raw)/length(meanFR120Same_Hab1Raw);
std_mean120Hab1_SameRaw = std(meanFR120Same_Hab1Raw);
n_mean120Hab1_SameRaw = length(meanFR120Same_Hab1Raw);

mean120Hab2_DiffRaw = sum(meanFR120Diff_Hab2Raw)/length(meanFR120Diff_Hab2Raw);
std_mean120Hab2_DiffRaw = std(meanFR120Diff_Hab2Raw);
n_mean120Hab2_DiffRaw = length(meanFR120Diff_Hab2Raw);
mean120Hab2_SameRaw = sum(meanFR120Same_Hab2Raw)/length(meanFR120Same_Hab2Raw);
std_mean120Hab2_SameRaw = std(meanFR120Same_Hab2Raw);
n_mean120Hab2_SameRaw = length(meanFR120Same_Hab2Raw);

mean120Trial1_DiffRaw = sum(meanFR120Diff_Trial1Raw)/length(meanFR120Diff_Trial1Raw);
std_mean120Trial1_DiffRaw = std(meanFR120Diff_Trial1Raw);
n_mean120Trial1_DiffRaw = length(meanFR120Diff_Trial1Raw);
mean120Trial1_SameRaw = sum(meanFR120Same_Trial1Raw)/length(meanFR120Same_Trial1Raw);
std_mean120Trial1_SameRaw = std(meanFR120Same_Trial1Raw);
n_mean120Trial1_SameRaw = length(meanFR120Same_Trial1Raw);

mean120Trial2_DiffRaw = sum(meanFR120Diff_Trial2Raw)/length(meanFR120Diff_Trial2Raw);
std_mean120Trial2_DiffRaw = std(meanFR120Diff_Trial2Raw);
n_mean120Trial2_DiffRaw = length(meanFR120Diff_Trial2Raw);
mean120Trial2_SameRaw = sum(meanFR120Same_Trial2Raw)/length(meanFR120Same_Trial2Raw);
std_mean120Trial2_SameRaw = std(meanFR120Same_Trial2Raw);
n_mean120Trial2_SameRaw = length(meanFR120Same_Trial2Raw);

meanFRstats120minsRaw = [mean120Hab1_DiffRaw,std_mean120Hab1_DiffRaw,n_mean120Hab1_DiffRaw;...
    mean120Hab1_SameRaw,std_mean120Hab1_SameRaw,n_mean120Hab1_SameRaw;...
    mean120Trial1_DiffRaw,std_mean120Trial1_DiffRaw,n_mean120Trial1_DiffRaw;...
    mean120Trial1_SameRaw,std_mean120Trial1_SameRaw,n_mean120Trial1_SameRaw;...
    mean120Hab2_DiffRaw,std_mean120Hab2_DiffRaw,n_mean120Hab2_DiffRaw;...
    mean120Hab2_SameRaw,std_mean120Hab2_SameRaw,n_mean120Hab2_SameRaw;...
    mean120Trial2_DiffRaw,std_mean120Trial2_DiffRaw,n_mean120Trial2_DiffRaw;...
    mean120Trial2_SameRaw,std_mean120Trial2_SameRaw,n_mean120Trial2_SameRaw;];

mean30Hab1_DiffRaw = sum(meanFR30Diff_Hab1Raw)/length(meanFR30Diff_Hab1Raw);
std_mean30Hab1_DiffRaw = std(meanFR30Diff_Hab1Raw);
n_mean30Hab1_DiffRaw = length(meanFR30Diff_Hab1Raw);
mean30Hab1_SameRaw = sum(meanFR30Same_Hab1Raw)/length(meanFR30Same_Hab1Raw);
std_mean30Hab1_SameRaw = std(meanFR30Same_Hab1Raw);
n_mean30Hab1_SameRaw = length(meanFR30Same_Hab1Raw);

mean30Hab2_DiffRaw = sum(meanFR30Diff_Hab2Raw)/length(meanFR30Diff_Hab2Raw);
std_mean30Hab2_DiffRaw = std(meanFR30Diff_Hab2Raw);
n_mean30Hab2_DiffRaw = length(meanFR30Diff_Hab2Raw);
mean30Hab2_SameRaw = sum(meanFR30Same_Hab2Raw)/length(meanFR30Same_Hab2Raw);
std_mean30Hab2_SameRaw = std(meanFR30Same_Hab2Raw);
n_mean30Hab2_SameRaw = length(meanFR30Same_Hab2Raw);

mean30Trial1_DiffRaw = sum(meanFR30Diff_Trial1Raw)/length(meanFR30Diff_Trial1Raw);
std_mean30Trial1_DiffRaw = std(meanFR30Diff_Trial1Raw);
n_mean30Trial1_DiffRaw = length(meanFR30Diff_Trial1Raw);
mean30Trial1_SameRaw = sum(meanFR30Same_Trial1Raw)/length(meanFR30Same_Trial1Raw);
std_mean30Trial1_SameRaw = std(meanFR30Same_Trial1Raw);
n_mean30Trial1_SameRaw = length(meanFR30Same_Trial1Raw);

mean30Trial2_DiffRaw = sum(meanFR30Diff_Trial2Raw)/length(meanFR30Diff_Trial2Raw);
std_mean30Trial2_DiffRaw = std(meanFR30Diff_Trial2Raw);
n_mean30Trial2_DiffRaw = length(meanFR30Diff_Trial2Raw);
mean30Trial2_SameRaw = sum(meanFR30Same_Trial2Raw)/length(meanFR30Same_Trial2Raw);
std_mean30Trial2_SameRaw = std(meanFR30Same_Trial2Raw);
n_mean30Trial2_SameRaw = length(meanFR30Same_Trial2Raw);

meanFRstats30minsRaw = [mean30Hab1_DiffRaw,std_mean30Hab1_DiffRaw,n_mean30Hab1_DiffRaw;...
    mean30Hab1_SameRaw,std_mean30Hab1_SameRaw,n_mean30Hab1_SameRaw;...
    mean30Trial1_DiffRaw,std_mean30Trial1_DiffRaw,n_mean30Trial1_DiffRaw;...
    mean30Trial1_SameRaw,std_mean30Trial1_SameRaw,n_mean30Trial1_SameRaw;...
    mean30Hab2_DiffRaw,std_mean30Hab2_DiffRaw,n_mean30Hab2_DiffRaw;...
    mean30Hab2_SameRaw,std_mean30Hab2_SameRaw,n_mean30Hab2_SameRaw;...
    mean30Trial2_DiffRaw,std_mean30Trial2_DiffRaw,n_mean30Trial2_DiffRaw;...
    mean30Trial2_SameRaw,std_mean30Trial2_SameRaw,n_mean30Trial2_SameRaw;];

% All Cells (no cell reg)

mean120Hab1_DiffAll = sum(meanFR120Diff_Hab1All)/length(meanFR120Diff_Hab1All);
std_mean120Hab1_DiffAll = std(meanFR120Diff_Hab1All);
n_mean120Hab1_DiffAll = length(meanFR120Diff_Hab1All);
mean120Hab1_SameAll = sum(meanFR120Same_Hab1All)/length(meanFR120Same_Hab1All);
std_mean120Hab1_SameAll = std(meanFR120Same_Hab1All);
n_mean120Hab1_SameAll = length(meanFR120Same_Hab1All);

mean120Hab2_DiffAll = sum(meanFR120Diff_Hab2All)/length(meanFR120Diff_Hab2All);
std_mean120Hab2_DiffAll = std(meanFR120Diff_Hab2All);
n_mean120Hab2_DiffAll = length(meanFR120Diff_Hab2All);
mean120Hab2_SameAll = sum(meanFR120Same_Hab2All)/length(meanFR120Same_Hab2All);
std_mean120Hab2_SameAll = std(meanFR120Same_Hab2All);
n_mean120Hab2_SameAll = length(meanFR120Same_Hab2All);

mean120Trial1_DiffAll = sum(meanFR120Diff_Trial1All)/length(meanFR120Diff_Trial1All);
std_mean120Trial1_DiffAll = std(meanFR120Diff_Trial1All);
n_mean120Trial1_DiffAll = length(meanFR120Diff_Trial1All);
mean120Trial1_SameAll = sum(meanFR120Same_Trial1All)/length(meanFR120Same_Trial1All);
std_mean120Trial1_SameAll = std(meanFR120Same_Trial1All);
n_mean120Trial1_SameAll = length(meanFR120Same_Trial1All);

mean120Trial2_DiffAll = sum(meanFR120Diff_Trial2All)/length(meanFR120Diff_Trial2All);
std_mean120Trial2_DiffAll = std(meanFR120Diff_Trial2All);
n_mean120Trial2_DiffAll = length(meanFR120Diff_Trial2All);
mean120Trial2_SameAll = sum(meanFR120Same_Trial2All)/length(meanFR120Same_Trial2All);
std_mean120Trial2_SameAll = std(meanFR120Same_Trial2All);
n_mean120Trial2_SameAll = length(meanFR120Same_Trial2All);

meanFRstats120minsAll = [mean120Hab1_DiffAll,std_mean120Hab1_DiffAll,n_mean120Hab1_DiffAll;...
    mean120Hab1_SameAll,std_mean120Hab1_SameAll,n_mean120Hab1_SameAll;...
    mean120Trial1_DiffAll,std_mean120Trial1_DiffAll,n_mean120Trial1_DiffAll;...
    mean120Trial1_SameAll,std_mean120Trial1_SameAll,n_mean120Trial1_SameAll;...
    mean120Hab2_DiffAll,std_mean120Hab2_DiffAll,n_mean120Hab2_DiffAll;...
    mean120Hab2_SameAll,std_mean120Hab2_SameAll,n_mean120Hab2_SameAll;...
    mean120Trial2_DiffAll,std_mean120Trial2_DiffAll,n_mean120Trial2_DiffAll;...
    mean120Trial2_SameAll,std_mean120Trial2_SameAll,n_mean120Trial2_SameAll;];

mean30Hab1_DiffAll = sum(meanFR30Diff_Hab1All)/length(meanFR30Diff_Hab1All);
std_mean30Hab1_DiffAll = std(meanFR30Diff_Hab1All);
n_mean30Hab1_DiffAll = length(meanFR30Diff_Hab1All);
mean30Hab1_SameAll = sum(meanFR30Same_Hab1All)/length(meanFR30Same_Hab1All);
std_mean30Hab1_SameAll = std(meanFR30Same_Hab1All);
n_mean30Hab1_SameAll = length(meanFR30Same_Hab1All);

mean30Hab2_DiffAll = sum(meanFR30Diff_Hab2All)/length(meanFR30Diff_Hab2All);
std_mean30Hab2_DiffAll = std(meanFR30Diff_Hab2All);
n_mean30Hab2_DiffAll = length(meanFR30Diff_Hab2All);
mean30Hab2_SameAll = sum(meanFR30Same_Hab2All)/length(meanFR30Same_Hab2All);
std_mean30Hab2_SameAll = std(meanFR30Same_Hab2All);
n_mean30Hab2_SameAll = length(meanFR30Same_Hab2All);

mean30Trial1_DiffAll = sum(meanFR30Diff_Trial1All)/length(meanFR30Diff_Trial1All);
std_mean30Trial1_DiffAll = std(meanFR30Diff_Trial1All);
n_mean30Trial1_DiffAll = length(meanFR30Diff_Trial1All);
mean30Trial1_SameAll = sum(meanFR30Same_Trial1All)/length(meanFR30Same_Trial1All);
std_mean30Trial1_SameAll = std(meanFR30Same_Trial1All);
n_mean30Trial1_SameAll = length(meanFR30Same_Trial1All);

mean30Trial2_DiffAll = sum(meanFR30Diff_Trial2All)/length(meanFR30Diff_Trial2All);
std_mean30Trial2_DiffAll = std(meanFR30Diff_Trial2All);
n_mean30Trial2_DiffAll = length(meanFR30Diff_Trial2All);
mean30Trial2_SameAll = sum(meanFR30Same_Trial2All)/length(meanFR30Same_Trial2All);
std_mean30Trial2_SameAll = std(meanFR30Same_Trial2All);
n_mean30Trial2_SameAll = length(meanFR30Same_Trial2All);

meanFRstats30minsAll = [mean30Hab1_DiffAll,std_mean30Hab1_DiffAll,n_mean30Hab1_DiffAll;...
    mean30Hab1_SameAll,std_mean30Hab1_SameAll,n_mean30Hab1_SameAll;...
    mean30Trial1_DiffAll,std_mean30Trial1_DiffAll,n_mean30Trial1_DiffAll;...
    mean30Trial1_SameAll,std_mean30Trial1_SameAll,n_mean30Trial1_SameAll;...
    mean30Hab2_DiffAll,std_mean30Hab2_DiffAll,n_mean30Hab2_DiffAll;...
    mean30Hab2_SameAll,std_mean30Hab2_SameAll,n_mean30Hab2_SameAll;...
    mean30Trial2_DiffAll,std_mean30Trial2_DiffAll,n_mean30Trial2_DiffAll;...
    mean30Trial2_SameAll,std_mean30Trial2_SameAll,n_mean30Trial2_SameAll;];

% passed only z nomralized

mean120Hab1_DiffAll_passed = sum(meanFR120Diff_Hab1All_passed)/length(meanFR120Diff_Hab1All_passed);
std_mean120Hab1_DiffAll_passed = std(meanFR120Diff_Hab1All_passed);
n_mean120Hab1_DiffAll_passed = length(meanFR120Diff_Hab1All_passed);
mean120Hab1_SameAll_passed = sum(meanFR120Same_Hab1All_passed)/length(meanFR120Same_Hab1All_passed);
std_mean120Hab1_SameAll_passed = std(meanFR120Same_Hab1All_passed);
n_mean120Hab1_SameAll_passed = length(meanFR120Same_Hab1All_passed);

mean120Hab2_DiffAll_passed = sum(meanFR120Diff_Hab2All_passed)/length(meanFR120Diff_Hab2All_passed);
std_mean120Hab2_DiffAll_passed = std(meanFR120Diff_Hab2All_passed);
n_mean120Hab2_DiffAll_passed = length(meanFR120Diff_Hab2All_passed);
mean120Hab2_SameAll_passed = sum(meanFR120Same_Hab2All_passed)/length(meanFR120Same_Hab2All_passed);
std_mean120Hab2_SameAll_passed = std(meanFR120Same_Hab2All_passed);
n_mean120Hab2_SameAll_passed = length(meanFR120Same_Hab2All_passed);

mean120Trial1_DiffAll_passed = sum(meanFR120Diff_Trial1All_passed)/length(meanFR120Diff_Trial1All_passed);
std_mean120Trial1_DiffAll_passed = std(meanFR120Diff_Trial1All_passed);
n_mean120Trial1_DiffAll_passed = length(meanFR120Diff_Trial1All_passed);
mean120Trial1_SameAll_passed = sum(meanFR120Same_Trial1All_passed)/length(meanFR120Same_Trial1All_passed);
std_mean120Trial1_SameAll_passed = std(meanFR120Same_Trial1All_passed);
n_mean120Trial1_SameAll_passed = length(meanFR120Same_Trial1All_passed);

mean120Trial2_DiffAll_passed = sum(meanFR120Diff_Trial2All_passed)/length(meanFR120Diff_Trial2All_passed);
std_mean120Trial2_DiffAll_passed = std(meanFR120Diff_Trial2All_passed);
n_mean120Trial2_DiffAll_passed = length(meanFR120Diff_Trial2All_passed);
mean120Trial2_SameAll_passed = sum(meanFR120Same_Trial2All_passed)/length(meanFR120Same_Trial2All_passed);
std_mean120Trial2_SameAll_passed = std(meanFR120Same_Trial2All_passed);
n_mean120Trial2_SameAll_passed = length(meanFR120Same_Trial2All_passed);

meanFRstats120minsAll_passed = [mean120Hab1_DiffAll_passed,std_mean120Hab1_DiffAll_passed,n_mean120Hab1_DiffAll_passed;...
    mean120Hab1_SameAll_passed,std_mean120Hab1_SameAll_passed,n_mean120Hab1_SameAll_passed;...
    mean120Trial1_DiffAll_passed,std_mean120Trial1_DiffAll_passed,n_mean120Trial1_DiffAll_passed;...
    mean120Trial1_SameAll_passed,std_mean120Trial1_SameAll_passed,n_mean120Trial1_SameAll_passed;...
    mean120Hab2_DiffAll_passed,std_mean120Hab2_DiffAll_passed,n_mean120Hab2_DiffAll_passed;...
    mean120Hab2_SameAll_passed,std_mean120Hab2_SameAll_passed,n_mean120Hab2_SameAll_passed;...
    mean120Trial2_DiffAll_passed,std_mean120Trial2_DiffAll_passed,n_mean120Trial2_DiffAll_passed;...
    mean120Trial2_SameAll_passed,std_mean120Trial2_SameAll_passed,n_mean120Trial2_SameAll_passed;];

mean30Hab1_DiffAll_passed = sum(meanFR30Diff_Hab1All_passed)/length(meanFR30Diff_Hab1All_passed);
std_mean30Hab1_DiffAll_passed = std(meanFR30Diff_Hab1All_passed);
n_mean30Hab1_DiffAll_passed = length(meanFR30Diff_Hab1All_passed);
mean30Hab1_SameAll_passed = sum(meanFR30Same_Hab1All_passed)/length(meanFR30Same_Hab1All_passed);
std_mean30Hab1_SameAll_passed = std(meanFR30Same_Hab1All_passed);
n_mean30Hab1_SameAll_passed = length(meanFR30Same_Hab1All_passed);

mean30Hab2_DiffAll_passed = sum(meanFR30Diff_Hab2All_passed)/length(meanFR30Diff_Hab2All_passed);
std_mean30Hab2_DiffAll_passed = std(meanFR30Diff_Hab2All_passed);
n_mean30Hab2_DiffAll_passed = length(meanFR30Diff_Hab2All_passed);
mean30Hab2_SameAll_passed = sum(meanFR30Same_Hab2All_passed)/length(meanFR30Same_Hab2All_passed);
std_mean30Hab2_SameAll_passed = std(meanFR30Same_Hab2All_passed);
n_mean30Hab2_SameAll_passed = length(meanFR30Same_Hab2All_passed);

mean30Trial1_DiffAll_passed = sum(meanFR30Diff_Trial1All)/length(meanFR30Diff_Trial1All_passed);
std_mean30Trial1_DiffAll_passed = std(meanFR30Diff_Trial1All_passed);
n_mean30Trial1_DiffAll_passed = length(meanFR30Diff_Trial1All_passed);
mean30Trial1_SameAll_passed = sum(meanFR30Same_Trial1All_passed)/length(meanFR30Same_Trial1All_passed);
std_mean30Trial1_SameAll_passed = std(meanFR30Same_Trial1All_passed);
n_mean30Trial1_SameAll_passed = length(meanFR30Same_Trial1All_passed);

mean30Trial2_DiffAll_passed = sum(meanFR30Diff_Trial2All_passed)/length(meanFR30Diff_Trial2All_passed);
std_mean30Trial2_DiffAll_passed = std(meanFR30Diff_Trial2All_passed);
n_mean30Trial2_DiffAll_passed = length(meanFR30Diff_Trial2All_passed);
mean30Trial2_SameAll_passed = sum(meanFR30Same_Trial2All_passed)/length(meanFR30Same_Trial2All_passed);
std_mean30Trial2_SameAll_passed = std(meanFR30Same_Trial2All_passed);
n_mean30Trial2_SameAll_passed = length(meanFR30Same_Trial2All_passed);

meanFRstats30minsAll_passed = [mean30Hab1_DiffAll_passed,std_mean30Hab1_DiffAll_passed,n_mean30Hab1_DiffAll_passed;...
    mean30Hab1_SameAll_passed,std_mean30Hab1_SameAll_passed,n_mean30Hab1_SameAll_passed;...
    mean30Trial1_DiffAll_passed,std_mean30Trial1_DiffAll_passed,n_mean30Trial1_DiffAll_passed;...
    mean30Trial1_SameAll_passed,std_mean30Trial1_SameAll_passed,n_mean30Trial1_SameAll_passed;...
    mean30Hab2_DiffAll_passed,std_mean30Hab2_DiffAll_passed,n_mean30Hab2_DiffAll_passed;...
    mean30Hab2_SameAll_passed,std_mean30Hab2_SameAll_passed,n_mean30Hab2_SameAll_passed;...
    mean30Trial2_DiffAll_passed,std_mean30Trial2_DiffAll_passed,n_mean30Trial2_DiffAll_passed;...
    mean30Trial2_SameAll_passed,std_mean30Trial2_SameAll_passed,n_mean30Trial2_SameAll_passed;];


% raw
mean120Hab1_DiffRawAll = sum(meanFR120Diff_Hab1Raw)/length(meanFR120Diff_Hab1Raw);
std_mean120Hab1_DiffRawAll = std(meanFR120Diff_Hab1Raw);
n_mean120Hab1_DiffRawAll = length(meanFR120Diff_Hab1Raw);
mean120Hab1_SameRawAll = sum(meanFR120Same_Hab1Raw)/length(meanFR120Same_Hab1Raw);
std_mean120Hab1_SameRawAll = std(meanFR120Same_Hab1Raw);
n_mean120Hab1_SameRawAll = length(meanFR120Same_Hab1Raw);

mean120Hab2_DiffRawAll = sum(meanFR120Diff_Hab2Raw)/length(meanFR120Diff_Hab2Raw);
std_mean120Hab2_DiffRawAll = std(meanFR120Diff_Hab2Raw);
n_mean120Hab2_DiffRawAll = length(meanFR120Diff_Hab2Raw);
mean120Hab2_SameRawAll = sum(meanFR120Same_Hab2Raw)/length(meanFR120Same_Hab2Raw);
std_mean120Hab2_SameRawAll = std(meanFR120Same_Hab2Raw);
n_mean120Hab2_SameRawAll = length(meanFR120Same_Hab2Raw);

mean120Trial1_DiffRawAll = sum(meanFR120Diff_Trial1Raw)/length(meanFR120Diff_Trial1Raw);
std_mean120Trial1_DiffRawAll = std(meanFR120Diff_Trial1Raw);
n_mean120Trial1_DiffRawAll = length(meanFR120Diff_Trial1Raw);
mean120Trial1_SameRawAll = sum(meanFR120Same_Trial1Raw)/length(meanFR120Same_Trial1Raw);
std_mean120Trial1_SameRawAll = std(meanFR120Same_Trial1Raw);
n_mean120Trial1_SameRawAll = length(meanFR120Same_Trial1Raw);

mean120Trial2_DiffRawAll = sum(meanFR120Diff_Trial2Raw)/length(meanFR120Diff_Trial2Raw);
std_mean120Trial2_DiffRawAll = std(meanFR120Diff_Trial2Raw);
n_mean120Trial2_DiffRawAll = length(meanFR120Diff_Trial2Raw);
mean120Trial2_SameRawAll = sum(meanFR120Same_Trial2Raw)/length(meanFR120Same_Trial2Raw);
std_mean120Trial2_SameRawAll = std(meanFR120Same_Trial2Raw);
n_mean120Trial2_SameRawAll = length(meanFR120Same_Trial2Raw);

meanFRstats120minsRawAll = [mean120Hab1_DiffRawAll,std_mean120Hab1_DiffRawAll,n_mean120Hab1_DiffRawAll;...
    mean120Hab1_SameRawAll,std_mean120Hab1_SameRawAll,n_mean120Hab1_SameRawAll;...
    mean120Trial1_DiffRawAll,std_mean120Trial1_DiffRawAll,n_mean120Trial1_DiffRawAll;...
    mean120Trial1_SameRawAll,std_mean120Trial1_SameRawAll,n_mean120Trial1_SameRawAll;...
    mean120Hab2_DiffRawAll,std_mean120Hab2_DiffRawAll,n_mean120Hab2_DiffRawAll;...
    mean120Hab2_SameRawAll,std_mean120Hab2_SameRawAll,n_mean120Hab2_SameRawAll;...
    mean120Trial2_DiffRawAll,std_mean120Trial2_DiffRawAll,n_mean120Trial2_DiffRawAll;...
    mean120Trial2_SameRawAll,std_mean120Trial2_SameRawAll,n_mean120Trial2_SameRawAll;];

mean30Hab1_DiffRawAll = sum(meanFR30Diff_Hab1Raw)/length(meanFR30Diff_Hab1Raw);
std_mean30Hab1_DiffRawAll = std(meanFR30Diff_Hab1Raw);
n_mean30Hab1_DiffRawAll = length(meanFR30Diff_Hab1Raw);
mean30Hab1_SameRawAll = sum(meanFR30Same_Hab1Raw)/length(meanFR30Same_Hab1Raw);
std_mean30Hab1_SameRawAll = std(meanFR30Same_Hab1Raw);
n_mean30Hab1_SameRawAll = length(meanFR30Same_Hab1Raw);

mean30Hab2_DiffRawAll = sum(meanFR30Diff_Hab2Raw)/length(meanFR30Diff_Hab2Raw);
std_mean30Hab2_DiffRawAll = std(meanFR30Diff_Hab2Raw);
n_mean30Hab2_DiffRawAll = length(meanFR30Diff_Hab2Raw);
mean30Hab2_SameRawAll = sum(meanFR30Same_Hab2Raw)/length(meanFR30Same_Hab2Raw);
std_mean30Hab2_SameRawAll = std(meanFR30Same_Hab2Raw);
n_mean30Hab2_SameRawAll = length(meanFR30Same_Hab2Raw);

mean30Trial1_DiffRawAll = sum(meanFR30Diff_Trial1Raw)/length(meanFR30Diff_Trial1Raw);
std_mean30Trial1_DiffRawAll = std(meanFR30Diff_Trial1Raw);
n_mean30Trial1_DiffRawAll = length(meanFR30Diff_Trial1Raw);
mean30Trial1_SameRawAll = sum(meanFR30Same_Trial1Raw)/length(meanFR30Same_Trial1Raw);
std_mean30Trial1_SameRawAll = std(meanFR30Same_Trial1Raw);
n_mean30Trial1_SameRawAll = length(meanFR30Same_Trial1Raw);

mean30Trial2_DiffRawAll = sum(meanFR30Diff_Trial2Raw)/length(meanFR30Diff_Trial2Raw);
std_mean30Trial2_DiffRawAll = std(meanFR30Diff_Trial2Raw);
n_mean30Trial2_DiffRawAll = length(meanFR30Diff_Trial2Raw);
mean30Trial2_SameRawAll = sum(meanFR30Same_Trial2Raw)/length(meanFR30Same_Trial2Raw);
std_mean30Trial2_SameRawAll = std(meanFR30Same_Trial2Raw);
n_mean30Trial2_SameRawAll = length(meanFR30Same_Trial2Raw);

meanFRstats30minsRawAll = [mean30Hab1_DiffRawAll,std_mean30Hab1_DiffRawAll,n_mean30Hab1_DiffRawAll;...
    mean30Hab1_SameRawAll,std_mean30Hab1_SameRawAll,n_mean30Hab1_SameRawAll;...
    mean30Trial1_DiffRawAll,std_mean30Trial1_DiffRawAll,n_mean30Trial1_DiffRawAll;...
    mean30Trial1_SameRawAll,std_mean30Trial1_SameRawAll,n_mean30Trial1_SameRawAll;...
    mean30Hab2_DiffRawAll,std_mean30Hab2_DiffRawAll,n_mean30Hab2_DiffRawAll;...
    mean30Hab2_SameRawAll,std_mean30Hab2_SameRawAll,n_mean30Hab2_SameRawAll;...
    mean30Trial2_DiffRawAll,std_mean30Trial2_DiffRawAll,n_mean30Trial2_DiffRawAll;...
    mean30Trial2_SameRawAll,std_mean30Trial2_SameRawAll,n_mean30Trial2_SameRawAll;];

%% figure 6b Before/After interaction

meanbout1_2sec_30Same = mean(bout1_2sec_30Same,1);
meanbout2_2sec_30Same = mean(bout2_2sec_30Same,1);
meanboutHab1_2sec_30Same = mean(boutHab1_2sec_30Same,1);
meanboutHab2_2sec_30Same = mean(boutHab2_2sec_30Same,1);

meanbout1_2sec_30Diff = mean(bout1_2sec_30Diff,1);
meanbout2_2sec_30Diff = mean(bout2_2sec_30Diff,1);
meanboutHab1_2sec_30Diff = mean(boutHab1_2sec_30Diff,1);
meanboutHab2_2sec_30Diff = mean(boutHab2_2sec_30Diff,1);

meanbout1_2sec_120Same = mean(bout1_2sec_120Same,1);
meanbout2_2sec_120Same = mean(bout2_2sec_120Same,1);
meanboutHab1_2sec_120Same = mean(boutHab1_2sec_120Same,1);
meanboutHab2_2sec_120Same = mean(boutHab2_2sec_120Same,1);

meanbout1_2sec_120Diff = mean(bout1_2sec_120Diff,1);
meanbout2_2sec_120Diff = mean(bout2_2sec_120Diff,1);
meanboutHab1_2sec_120Diff = mean(boutHab1_2sec_120Diff,1);
meanboutHab2_2sec_120Diff = mean(boutHab2_2sec_120Diff,1);

meanbout1_3sec_30Same = mean(bout1_3sec_30Same,1);
meanbout2_3sec_30Same = mean(bout2_3sec_30Same,1);
meanboutHab1_3sec_30Same = mean(boutHab1_3sec_30Same,1);
meanboutHab2_3sec_30Same = mean(boutHab2_3sec_30Same,1);

meanbout1_3sec_30Diff = mean(bout1_3sec_30Diff,1);
meanbout2_3sec_30Diff = mean(bout2_3sec_30Diff,1);
meanboutHab1_3sec_30Diff = mean(boutHab1_3sec_30Diff,1);
meanboutHab2_3sec_30Diff = mean(boutHab2_3sec_30Diff,1);

meanbout1_3sec_120Same = mean(bout1_3sec_120Same,1);
meanbout2_3sec_120Same = mean(bout2_3sec_120Same,1);
meanboutHab1_3sec_120Same = mean(boutHab1_3sec_120Same,1);
meanboutHab2_3sec_120Same = mean(boutHab2_3sec_120Same,1);

meanbout1_3sec_120Diff = mean(bout1_3sec_120Diff,1);
meanbout2_3sec_120Diff = mean(bout2_3sec_120Diff,1);
meanboutHab1_3sec_120Diff = mean(boutHab1_3sec_120Diff,1);
meanboutHab2_3sec_120Diff = mean(boutHab2_3sec_120Diff,1);

meanbout1_32sec_120Diff = mean(bout1_32sec_120Diff,1);
meanbout2_32sec_120Diff = mean(bout2_32sec_120Diff,1);
meanboutHab1_32sec_120Diff = mean(boutHab1_32sec_120Diff,1);
meanboutHab2_32sec_120Diff = mean(boutHab2_32sec_120Diff,1);

meanbout1_32sec_120Same = mean(bout1_32sec_120Same,1);
meanbout2_32sec_120Same = mean(bout2_32sec_120Same,1);
meanboutHab1_32sec_120Same = mean(boutHab1_32sec_120Same,1);
meanboutHab2_32sec_120Same = mean(boutHab2_32sec_120Same,1);

meanbout1_32sec_30Diff = mean(bout1_32sec_30Diff,1);
meanbout2_32sec_30Diff = mean(bout2_32sec_30Diff,1);
meanboutHab1_32sec_30Diff = mean(boutHab1_32sec_30Diff,1);
meanboutHab2_32sec_30Diff = mean(boutHab2_32sec_30Diff,1);

meanbout1_32sec_30Same = mean(bout1_32sec_30Same,1);
meanbout2_32sec_30Same = mean(bout2_32sec_30Same,1);
meanboutHab1_32sec_30Same = mean(boutHab1_32sec_30Same,1);
meanboutHab2_32sec_30Same = mean(boutHab2_32sec_30Same,1);


meanboutHab1_2sec_30Same_1h = mean(boutHab1_2sec_30Same_1h,1);
meanboutHab2_2sec_30Same_1h = mean(boutHab2_2sec_30Same_1h,1);
meanboutHab1_2sec_30Diff_1h = mean(boutHab1_2sec_30Diff_1h,1);
meanboutHab2_2sec_30Diff_1h = mean(boutHab2_2sec_30Diff_1h,1);

meanboutHab1_2sec_120Same_1h = mean(boutHab1_2sec_120Same_1h,1);
meanboutHab2_2sec_120Same_1h = mean(boutHab2_2sec_120Same_1h,1);
meanboutHab1_2sec_120Diff_1h = mean(boutHab1_2sec_120Diff_1h,1);
meanboutHab2_2sec_120Diff_1h = mean(boutHab2_2sec_120Diff_1h,1);

meanboutHab1_2sec_30Same_2h = mean(boutHab1_2sec_30Same_2h,1);
meanboutHab2_2sec_30Same_2h = mean(boutHab2_2sec_30Same_2h,1);
meanboutHab1_2sec_30Diff_2h = mean(boutHab1_2sec_30Diff_2h,1);
meanboutHab2_2sec_30Diff_2h = mean(boutHab2_2sec_30Diff_2h,1);

meanboutHab1_2sec_120Same_2h = mean(boutHab1_2sec_120Same_2h,1);
meanboutHab2_2sec_120Same_2h = mean(boutHab2_2sec_120Same_2h,1);
meanboutHab1_2sec_120Diff_2h = mean(boutHab1_2sec_120Diff_2h,1);
meanboutHab2_2sec_120Diff_2h = mean(boutHab2_2sec_120Diff_2h,1);


[~,sortbout1_2sec_30Same] = sort(mean(bout1_2sec_30Same,2),'descend');
[~,sortbout2_2sec_30Same] = sort(mean(bout2_2sec_30Same,2),'descend');
[~,sortboutHab1_2sec_30Same] = sort(mean(boutHab1_2sec_30Same,2),'descend');
[~,sortboutHab2_2sec_30Same] = sort(mean(boutHab2_2sec_30Same,2),'descend');
[~,sortbout1_3sec_30Same] = sort(mean(bout1_3sec_30Same,2),'descend');
[~,sortbout2_3sec_30Same] = sort(mean(bout2_3sec_30Same,2),'descend');
[~,sortboutHab1_3sec_30Same] = sort(mean(boutHab1_3sec_30Same,2),'descend');
[~,sortboutHab2_3sec_30Same] = sort(mean(boutHab2_3sec_30Same,2),'descend');
[~,sortbout1_2sec_120Same] = sort(mean(bout1_2sec_120Same,2),'descend');
[~,sortbout2_2sec_120Same] = sort(mean(bout2_2sec_120Same,2),'descend');
[~,sortboutHab1_2sec_120Same] = sort(mean(boutHab1_2sec_120Same,2),'descend');
[~,sortboutHab2_2sec_120Same] = sort(mean(boutHab2_2sec_120Same,2),'descend');
[~,sortbout1_3sec_120Same] = sort(mean(bout1_3sec_120Same,2),'descend');
[~,sortbout2_3sec_120Same] = sort(mean(bout2_3sec_120Same,2),'descend');
[~,sortboutHab1_3sec_120Same] = sort(mean(boutHab1_3sec_120Same,2),'descend');
[~,sortboutHab2_3sec_120Same] = sort(mean(boutHab2_3sec_120Same,2),'descend');
[~,sortbout1_2sec_30Diff] = sort(mean(bout1_2sec_30Diff,2),'descend');
[~,sortbout2_2sec_30Diff] = sort(mean(bout2_2sec_30Diff,2),'descend');
[~,sortboutHab1_2sec_30Diff] = sort(mean(boutHab1_2sec_30Diff,2),'descend');
[~,sortboutHab2_2sec_30Diff] = sort(mean(boutHab2_2sec_30Diff,2),'descend');
[~,sortbout1_3sec_30Diff] = sort(mean(bout1_3sec_30Diff,2),'descend');
[~,sortbout2_3sec_30Diff] = sort(mean(bout2_3sec_30Diff,2),'descend');
[~,sortboutHab1_3sec_30Diff] = sort(mean(boutHab1_3sec_30Diff,2),'descend');
[~,sortboutHab2_3sec_30Diff] = sort(mean(boutHab2_3sec_30Diff,2),'descend');
[~,sortbout1_2sec_120Diff] = sort(mean(bout1_2sec_120Diff,2),'descend');
[~,sortbout2_2sec_120Diff] = sort(mean(bout2_2sec_120Diff,2),'descend');
[~,sortboutHab1_2sec_120Diff] = sort(mean(boutHab1_2sec_120Diff,2),'descend');
[~,sortboutHab2_2sec_120Diff] = sort(mean(boutHab2_2sec_120Diff,2),'descend');
[~,sortbout1_3sec_120Diff] = sort(mean(bout1_3sec_120Diff,2),'descend');
[~,sortbout2_3sec_120Diff] = sort(mean(bout2_3sec_120Diff,2),'descend');
[~,sortboutHab1_3sec_120Diff] = sort(mean(boutHab1_3sec_120Diff,2),'descend');
[~,sortboutHab2_3sec_120Diff] = sort(mean(boutHab2_3sec_120Diff,2),'descend');

[~,sortbout1_32sec_30Same] = sort(mean(bout1_32sec_30Same,2),'descend');
[~,sortbout2_32sec_30Same] = sort(mean(bout2_32sec_30Same,2),'descend');
[~,sortboutHab1_32sec_30Same] = sort(mean(boutHab1_32sec_30Same,2),'descend');
[~,sortboutHab2_32sec_30Same] = sort(mean(boutHab2_32sec_30Same,2),'descend');
[~,sortbout1_32sec_30Diff] = sort(mean(bout1_32sec_30Diff,2),'descend');
[~,sortbout2_32sec_30Diff] = sort(mean(bout2_32sec_30Diff,2),'descend');
[~,sortboutHab1_32sec_30Diff] = sort(mean(boutHab1_32sec_30Diff,2),'descend');
[~,sortboutHab2_32sec_30Diff] = sort(mean(boutHab2_32sec_120Diff,2),'descend');
[~,sortbout1_32sec_120Same] = sort(mean(bout1_32sec_120Same,2),'descend');
[~,sortbout2_32sec_120Same] = sort(mean(bout2_32sec_120Same,2),'descend');
[~,sortboutHab1_32sec_120Same] = sort(mean(boutHab1_32sec_120Same,2),'descend');
[~,sortboutHab2_32sec_120Same] = sort(mean(boutHab2_32sec_120Same,2),'descend');
[~,sortbout1_32sec_120Diff] = sort(mean(bout1_32sec_120Diff,2),'descend');
[~,sortbout2_32sec_120Diff] = sort(mean(bout2_32sec_120Diff,2),'descend');
[~,sortboutHab1_32sec_120Diff] = sort(mean(boutHab1_32sec_120Diff,2),'descend');
[~,sortboutHab2_32sec_120Diff] = sort(mean(boutHab2_32sec_120Diff,2),'descend');

[~,sortboutHab1_2sec_30Same_1h] = sort(mean(boutHab1_2sec_30Same_1h,2),'descend');
[~,sortboutHab2_2sec_30Same_1h] = sort(mean(boutHab2_2sec_30Same_1h,2),'descend');
[~,sortboutHab1_2sec_120Same_1h] = sort(mean(boutHab1_2sec_120Same_1h,2),'descend');
[~,sortboutHab2_2sec_120Same_1h] = sort(mean(boutHab2_2sec_120Same_1h,2),'descend');
[~,sortboutHab1_2sec_30Diff_1h] = sort(mean(boutHab1_2sec_30Diff_1h,2),'descend');
[~,sortboutHab2_2sec_30Diff_1h] = sort(mean(boutHab2_2sec_30Diff_1h,2),'descend');
[~,sortboutHab1_2sec_120Diff_1h] = sort(mean(boutHab1_2sec_120Diff_1h,2),'descend');
[~,sortboutHab2_2sec_120Diff_1h] = sort(mean(boutHab2_2sec_120Diff_1h,2),'descend');

[~,sortboutHab1_2sec_30Same_2h] = sort(mean(boutHab1_2sec_30Same_2h,2),'descend');
[~,sortboutHab2_2sec_30Same_2h] = sort(mean(boutHab2_2sec_30Same_2h,2),'descend');
[~,sortboutHab1_2sec_120Same_2h] = sort(mean(boutHab1_2sec_120Same_2h,2),'descend');
[~,sortboutHab2_2sec_120Same_2h] = sort(mean(boutHab2_2sec_120Same_2h,2),'descend');
[~,sortboutHab1_2sec_30Diff_2h] = sort(mean(boutHab1_2sec_30Diff_2h,2),'descend');
[~,sortboutHab2_2sec_30Diff_2h] = sort(mean(boutHab2_2sec_30Diff_2h,2),'descend');
[~,sortboutHab1_2sec_120Diff_2h] = sort(mean(boutHab1_2sec_120Diff_2h,2),'descend');
[~,sortboutHab2_2sec_120Diff_2h] = sort(mean(boutHab2_2sec_120Diff_2h,2),'descend');

minbout2sec = min(min(cat(1,bout1_2sec_30Same,bout1_2sec_30Diff,bout1_2sec_120Same,bout1_2sec_120Diff,bout2_2sec_30Same,bout2_2sec_30Diff,bout2_2sec_120Same,bout2_2sec_120Diff,boutHab1_2sec_30Same,boutHab1_2sec_30Diff,boutHab1_2sec_120Same,boutHab1_2sec_120Diff,boutHab2_2sec_30Same,boutHab2_2sec_30Diff,boutHab2_2sec_120Same,boutHab2_2sec_120Diff)));
maxbout2sec = max(max(cat(1,bout1_2sec_30Same,bout1_2sec_30Diff,bout1_2sec_120Same,bout1_2sec_120Diff,bout2_2sec_30Same,bout2_2sec_30Diff,bout2_2sec_120Same,bout2_2sec_120Diff,boutHab1_2sec_30Same,boutHab1_2sec_30Diff,boutHab1_2sec_120Same,boutHab1_2sec_120Diff,boutHab2_2sec_30Same,boutHab2_2sec_30Diff,boutHab2_2sec_120Same,boutHab2_2sec_120Diff)));

minbout3sec = min(min(cat(1,bout1_3sec_30Same,bout1_3sec_30Diff,bout1_3sec_120Same,bout1_3sec_120Diff,bout2_3sec_30Same,bout2_3sec_30Diff,bout2_3sec_120Same,bout2_3sec_120Diff,boutHab1_3sec_30Same,boutHab1_3sec_30Diff,boutHab1_3sec_120Same,boutHab1_3sec_120Diff,boutHab2_3sec_30Same,boutHab2_3sec_30Diff,boutHab2_3sec_120Same,boutHab2_3sec_120Diff)));
maxbout3sec = max(max(cat(1,bout1_3sec_30Same,bout1_3sec_30Diff,bout1_3sec_120Same,bout1_3sec_120Diff,bout2_3sec_30Same,bout2_3sec_30Diff,bout2_3sec_120Same,bout2_3sec_120Diff,boutHab1_3sec_30Same,boutHab1_3sec_30Diff,boutHab1_3sec_120Same,boutHab1_3sec_120Diff,boutHab2_3sec_30Same,boutHab2_3sec_30Diff,boutHab2_3sec_120Same,boutHab2_3sec_120Diff)));

minboutmean2sec = min(min(cat(1,meanbout1_2sec_30Same,meanbout1_2sec_30Diff,meanbout1_2sec_120Same,meanbout1_2sec_120Diff,meanbout2_2sec_30Same,meanbout2_2sec_30Diff,meanbout2_2sec_120Same,meanbout2_2sec_120Diff,meanboutHab1_2sec_30Same,meanboutHab1_2sec_30Diff,meanboutHab1_2sec_120Same,meanboutHab1_2sec_120Diff,meanboutHab2_2sec_30Same,meanboutHab2_2sec_30Diff,meanboutHab2_2sec_120Same,meanboutHab2_2sec_120Diff)));
maxboutmean2sec = max(max(cat(1,meanbout1_2sec_30Same,meanbout1_2sec_30Diff,meanbout1_2sec_120Same,meanbout1_2sec_120Diff,meanbout2_2sec_30Same,meanbout2_2sec_30Diff,meanbout2_2sec_120Same,meanbout2_2sec_120Diff,meanboutHab1_2sec_30Same,meanboutHab1_2sec_30Diff,meanboutHab1_2sec_120Same,meanboutHab1_2sec_120Diff,meanboutHab2_2sec_30Same,meanboutHab2_2sec_30Diff,meanboutHab2_2sec_120Same,meanboutHab2_2sec_120Diff)));

minboutmean3sec = min(min(cat(1,meanbout1_3sec_30Same,meanbout1_3sec_30Diff,meanbout1_3sec_120Same,meanbout1_3sec_120Diff,meanbout2_3sec_30Same,meanbout2_3sec_30Diff,meanbout2_3sec_120Same,meanbout2_3sec_120Diff,meanboutHab1_3sec_30Same,meanboutHab1_3sec_30Diff,meanboutHab1_3sec_120Same,meanboutHab1_3sec_120Diff,meanboutHab2_3sec_30Same,meanboutHab2_3sec_30Diff,meanboutHab2_3sec_120Same,meanboutHab2_3sec_120Diff)));
maxboutmean3sec = max(max(cat(1,meanbout1_3sec_30Same,meanbout1_3sec_30Diff,meanbout1_3sec_120Same,meanbout1_3sec_120Diff,meanbout2_3sec_30Same,meanbout2_3sec_30Diff,meanbout2_3sec_120Same,meanbout2_3sec_120Diff,meanboutHab1_3sec_30Same,meanboutHab1_3sec_30Diff,meanboutHab1_3sec_120Same,meanboutHab1_3sec_120Diff,meanboutHab2_3sec_30Same,meanboutHab2_3sec_30Diff,meanboutHab2_3sec_120Same,meanboutHab2_3sec_120Diff)));

%30 mins Same 2 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_2sec_30Same(sortbout1_2sec_30Same,:))
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([minbout2sec maxbout2sec])
colorbar
title('Trial1 2-seconds')
subplot(2,1,2)
plot(meanbout1_2sec_30Same)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_30Same_2sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_2sec_30Same(sortbout2_2sec_30Same,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Trial2 2-seconds')
subplot(2,1,2)
plot(meanbout2_2sec_30Same)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_30Same_2sec')
% saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Deconvolved.jpeg')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_30Same(sortboutHab1_2sec_30Same,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds')
subplot(2,1,2)
plot(meanboutHab1_2sec_30Same)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Same_2sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_30Same(sortboutHab2_2sec_30Same,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds')
subplot(2,1,2)
plot(meanboutHab2_2sec_30Same)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Same_2sec')

%30 mins Same 3 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_3sec_30Same(sortbout1_3sec_30Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial1 3-seconds')
subplot(2,1,2)
plot(meanbout1_3sec_30Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Familiar')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Trial1_30Same_3sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_3sec_30Same(sortbout2_3sec_30Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial2 3-seconds')
subplot(2,1,2)
plot(meanbout2_3sec_30Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Familiar')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Trial2_30Same_3sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_3sec_30Same(sortboutHab1_3sec_30Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 1 3-seconds')
subplot(2,1,2)
plot(meanboutHab1_3sec_30Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Familiar')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Same_3sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_3sec_30Same(sortboutHab2_3sec_30Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 2 3-seconds')
subplot(2,1,2)
plot(meanboutHab2_3sec_30Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Familiar')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Same_3sec')

%%%%%%%%%%%%%
%120 mins Same 2 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_2sec_120Same(sortbout1_2sec_120Same,:))
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([minbout2sec maxbout2sec])
colorbar
title('Trial1 2-seconds')
subplot(2,1,2)
plot(meanbout1_2sec_120Same)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_120Same_2sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_2sec_120Same(sortbout2_2sec_120Same,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Trial2 2-seconds')
subplot(2,1,2)
plot(meanbout2_2sec_120Same)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_120Same_2sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_120Same(sortboutHab1_2sec_120Same,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds')
subplot(2,1,2)
plot(meanboutHab1_2sec_120Same)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Same_2sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_120Same(sortboutHab2_2sec_120Same,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds')
subplot(2,1,2)
plot(meanboutHab2_2sec_120Same)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Same_2sec')

%120 mins Same 3 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_3sec_120Same(sortbout1_3sec_120Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial1 3-seconds')
subplot(2,1,2)
plot(meanbout1_3sec_120Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Familiar')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Trial1_120Same_3sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_3sec_120Same(sortbout2_3sec_120Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial2 3-seconds')
subplot(2,1,2)
plot(meanbout2_3sec_120Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Familiar')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Trial2_120Same_3sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_3sec_120Same(sortboutHab1_3sec_120Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 1 3-seconds')
subplot(2,1,2)
plot(meanboutHab1_3sec_120Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Familiar')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Same_3sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_3sec_120Same(sortboutHab2_3sec_120Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 2 3-seconds')
subplot(2,1,2)
plot(meanboutHab2_3sec_120Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Familiar')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Same_3sec')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%30 mins Diff 2 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_2sec_30Diff(sortbout1_2sec_30Diff,:))
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([minbout2sec maxbout2sec])
colorbar
title('Trial1 2-seconds')
subplot(2,1,2)
plot(meanbout1_2sec_30Diff)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_30Diff_2sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_2sec_30Diff(sortbout2_2sec_30Diff,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Trial2 2-seconds')
subplot(2,1,2)
plot(meanbout2_2sec_30Diff)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_30Diff_2sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_30Diff(sortboutHab1_2sec_30Diff,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds')
subplot(2,1,2)
plot(meanboutHab1_2sec_30Diff)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Diff_2sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_30Diff(sortboutHab2_2sec_30Diff,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds')
subplot(2,1,2)
plot(meanboutHab2_2sec_30Diff)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Diff_2sec')

%30 mins Same 3 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_3sec_30Diff(sortbout1_3sec_30Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial1 3-seconds')
subplot(2,1,2)
plot(meanbout1_3sec_30Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Novel')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Trial1_30Diff_3sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_3sec_30Diff(sortbout2_3sec_30Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial2 3-seconds')
subplot(2,1,2)
plot(meanbout2_3sec_30Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Novel')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Trial2_30Diff_3sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_3sec_30Diff(sortboutHab1_3sec_30Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 1 3-seconds')
subplot(2,1,2)
plot(meanboutHab1_3sec_30Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Novel')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Diff_3sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_3sec_30Diff(sortboutHab2_3sec_30Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 2 3-seconds')
subplot(2,1,2)
plot(meanboutHab2_3sec_30Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Novel')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Diff_3sec')

%%%%%%%%%%%%%
%120 mins Diff 2 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_2sec_120Diff(sortbout1_2sec_120Diff,:))
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([minbout2sec maxbout2sec])
colorbar
title('Trial1 2-seconds')
subplot(2,1,2)
plot(meanbout1_2sec_120Diff)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_120Diff_2sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_2sec_120Diff(sortbout2_2sec_120Diff,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Trial2 2-seconds')
subplot(2,1,2)
plot(meanbout2_2sec_120Diff)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_120Diff_2sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_120Diff(sortboutHab1_2sec_120Diff,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds')
subplot(2,1,2)
plot(meanboutHab1_2sec_120Diff)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Diff_2sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_120Diff(sortboutHab2_2sec_120Diff,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds')
subplot(2,1,2)
plot(meanboutHab2_2sec_120Diff)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Diff_2sec')

%120 mins Diff 3 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_3sec_120Diff(sortbout1_3sec_120Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial1 3-seconds')
subplot(2,1,2)
plot(meanbout1_3sec_120Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Novel')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Trial1_120Diff_3sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_3sec_120Diff(sortbout2_3sec_120Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial2 3-seconds')
subplot(2,1,2)
plot(meanbout2_3sec_120Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Novel')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Trial2_120Diff_3sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_3sec_120Diff(sortboutHab1_3sec_120Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 1 3-seconds')
subplot(2,1,2)
plot(meanboutHab1_3sec_120Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Novel')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Diff_3sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_3sec_120Diff(sortboutHab2_3sec_120Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1])
xticklabels({'-3','-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 2 3-seconds')
subplot(2,1,2)
plot(meanboutHab2_3sec_120Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Novel')
vline(secframes3,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Diff_3sec')

%120 mins Diff 2 pre 3 seconds post interaction
%trial1
figure
subplot(2,1,1)
imagesc(bout1_32sec_120Diff(sortbout1_32sec_120Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial1 3-seconds')
subplot(2,1,2)
plot(meanbout1_32sec_120Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_120Diff_3-2sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_3sec_120Diff(sortbout2_3sec_120Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial2 3-seconds')
subplot(2,1,2)
plot(meanbout2_32sec_120Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_120Diff_3-2sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_3sec_120Diff(sortboutHab1_3sec_120Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 1 3-seconds')
subplot(2,1,2)
plot(meanboutHab1_32sec_120Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Diff_3-2sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_3sec_120Diff(sortboutHab2_3sec_120Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 2 3-seconds')
subplot(2,1,2)
plot(meanboutHab2_32sec_120Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Diff_3-2sec')

% 30 mins Diff 2 pre 3 seconds post interaction
%trial1
figure
subplot(2,1,1)
imagesc(bout1_32sec_30Diff(sortbout1_32sec_30Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial1 3-seconds')
subplot(2,1,2)
plot(meanbout1_32sec_30Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_30Diff_3-2sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_3sec_30Diff(sortbout2_3sec_30Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial2 3-seconds')
subplot(2,1,2)
plot(meanbout2_32sec_30Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_30Diff_3-2sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_3sec_30Diff(sortboutHab1_3sec_30Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 1 3-seconds')
subplot(2,1,2)
plot(meanboutHab1_32sec_30Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Diff_3-2sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_3sec_30Diff(sortboutHab2_3sec_30Diff,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 2 3-seconds')
subplot(2,1,2)
plot(meanboutHab2_32sec_30Diff)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Diff_3-2sec')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%120 mins Same 2 pre 3 seconds post interaction
%trial1
figure
subplot(2,1,1)
imagesc(bout1_32sec_120Same(sortbout1_32sec_120Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial1 3-seconds')
subplot(2,1,2)
plot(meanbout1_32sec_120Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_120Same_3-2sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_3sec_120Same(sortbout2_3sec_120Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial2 3-seconds')
subplot(2,1,2)
plot(meanbout2_32sec_120Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_120Same_3-2sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_3sec_120Same(sortboutHab1_3sec_120Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 1 3-seconds')
subplot(2,1,2)
plot(meanboutHab1_32sec_120Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Same_3-2sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_3sec_120Same(sortboutHab2_3sec_120Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 2 3-seconds')
subplot(2,1,2)
plot(meanboutHab2_32sec_120Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Same_3-2sec')

% 30 mins Diff 2 pre 3 seconds post interaction
%trial1
figure
subplot(2,1,1)
imagesc(bout1_32sec_30Same(sortbout1_32sec_30Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial1 3-seconds')
subplot(2,1,2)
plot(meanbout1_32sec_30Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_30Same_3-2sec')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_3sec_30Same(sortbout2_3sec_30Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Trial2 3-seconds')
subplot(2,1,2)
plot(meanbout2_32sec_30Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_30Same_3-2sec')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_3sec_30Same(sortboutHab1_3sec_30Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 1 3-seconds')
subplot(2,1,2)
plot(meanboutHab1_32sec_30Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Same_3-2sec')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_3sec_30Same(sortboutHab2_3sec_30Same,:))
caxis([minbout3sec maxbout3sec])
xticks([1:30:secframes3*2+1-30])
xticklabels({'-2', '-1', '0', '1', '2','3'})
colorbar
title('Habituation 2 3-seconds')
subplot(2,1,2)
plot(meanboutHab2_32sec_30Same)
ylim([minboutmean3sec maxboutmean3sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Same_3-2sec')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure6 Lopez Cell Reg

meanbout1_2sec_30Same_CellReg = mean(bout1_2sec_30Same_CellReg,1);
meanbout2_2sec_30Same_CellReg = mean(bout2_2sec_30Same_CellReg,1);
meanboutHab1_2sec_30Same_CellReg = mean(boutHab1_2sec_30Same_CellReg,1);
meanboutHab2_2sec_30Same_CellReg = mean(boutHab2_2sec_30Same_CellReg,1);

meanbout1_2sec_30Diff_CellReg = mean(bout1_2sec_30Diff_CellReg,1);
meanbout2_2sec_30Diff_CellReg = mean(bout2_2sec_30Diff_CellReg,1);
meanboutHab1_2sec_30Diff_CellReg = mean(boutHab1_2sec_30Diff_CellReg,1);
meanboutHab2_2sec_30Diff_CellReg = mean(boutHab2_2sec_30Diff_CellReg,1);

meanbout1_2sec_120Same_CellReg = mean(bout1_2sec_120Same_CellReg,1);
meanbout2_2sec_120Same_CellReg = mean(bout2_2sec_120Same_CellReg,1);
meanboutHab1_2sec_120Same_CellReg = mean(boutHab1_2sec_120Same_CellReg,1);
meanboutHab2_2sec_120Same_CellReg = mean(boutHab2_2sec_120Same_CellReg,1);

meanbout1_2sec_120Diff_CellReg = mean(bout1_2sec_120Diff_CellReg,1);
meanbout2_2sec_120Diff_CellReg = mean(bout2_2sec_120Diff_CellReg,1);
meanboutHab1_2sec_120Diff_CellReg = mean(boutHab1_2sec_120Diff_CellReg,1);
meanboutHab2_2sec_120Diff_CellReg = mean(boutHab2_2sec_120Diff_CellReg,1);

[~,sortbout1_2sec_120Diff_CellReg] = sort(mean(bout1_2sec_120Diff_CellReg,2),'descend');
[~,sortbout2_2sec_120Diff_CellReg] = sort(mean(bout2_2sec_120Diff_CellReg,2),'descend');
[~,sortboutHab1_2sec_120Diff_CellReg] = sort(mean(boutHab1_2sec_120Diff_CellReg,2),'descend');
[~,sortboutHab2_2sec_120Diff_CellReg] = sort(mean(boutHab2_2sec_120Diff_CellReg,2),'descend');
[~,sortbout1_2sec_30Diff_CellReg] = sort(mean(bout1_2sec_30Diff_CellReg,2),'descend');
[~,sortbout2_2sec_30Diff_CellReg] = sort(mean(bout2_2sec_30Diff_CellReg,2),'descend');
[~,sortboutHab1_2sec_30Diff_CellReg] = sort(mean(boutHab1_2sec_30Diff_CellReg,2),'descend');
[~,sortboutHab2_2sec_30Diff_CellReg] = sort(mean(boutHab2_2sec_30Diff_CellReg,2),'descend');
[~,sortbout1_2sec_120Same_CellReg] = sort(mean(bout1_2sec_120Same_CellReg,2),'descend');
[~,sortbout2_2sec_120Same_CellReg] = sort(mean(bout2_2sec_120Same_CellReg,2),'descend');
[~,sortboutHab1_2sec_120Same_CellReg] = sort(mean(boutHab1_2sec_120Same_CellReg,2),'descend');
[~,sortboutHab2_2sec_120Same_CellReg] = sort(mean(boutHab2_2sec_120Same_CellReg,2),'descend');
[~,sortbout1_2sec_30Same_CellReg] = sort(mean(bout1_2sec_30Same_CellReg,2),'descend');
[~,sortbout2_2sec_30Same_CellReg] = sort(mean(bout2_2sec_30Same_CellReg,2),'descend');
[~,sortboutHab1_2sec_30Same_CellReg] = sort(mean(boutHab1_2sec_30Same_CellReg,2),'descend');
[~,sortboutHab2_2sec_30Same_CellReg] = sort(mean(boutHab2_2sec_30Same_CellReg,2),'descend');

minbout2sec_CellReg = min(min(cat(1,bout1_2sec_30Same_CellReg,bout1_2sec_30Diff_CellReg,bout1_2sec_120Same_CellReg,bout1_2sec_120Diff_CellReg,bout2_2sec_30Same_CellReg,bout2_2sec_30Diff_CellReg,bout2_2sec_120Same_CellReg,bout2_2sec_120Diff_CellReg,boutHab1_2sec_30Same_CellReg,boutHab1_2sec_30Diff_CellReg,boutHab1_2sec_120Same_CellReg,boutHab1_2sec_120Diff_CellReg,boutHab2_2sec_30Same_CellReg,boutHab2_2sec_30Diff_CellReg,boutHab2_2sec_120Same_CellReg,boutHab2_2sec_120Diff_CellReg)));
maxbout2sec_CellReg = max(max(cat(1,bout1_2sec_30Same_CellReg,bout1_2sec_30Diff_CellReg,bout1_2sec_120Same_CellReg,bout1_2sec_120Diff_CellReg,bout2_2sec_30Same_CellReg,bout2_2sec_30Diff_CellReg,bout2_2sec_120Same_CellReg,bout2_2sec_120Diff_CellReg,boutHab1_2sec_30Same_CellReg,boutHab1_2sec_30Diff_CellReg,boutHab1_2sec_120Same_CellReg,boutHab1_2sec_120Diff_CellReg,boutHab2_2sec_30Same_CellReg,boutHab2_2sec_30Diff_CellReg,boutHab2_2sec_120Same_CellReg,boutHab2_2sec_120Diff_CellReg)));

minboutmean2sec_CellReg = min(min(cat(1,meanbout1_2sec_30Same_CellReg,meanbout1_2sec_30Diff_CellReg,meanbout1_2sec_120Same_CellReg,meanbout1_2sec_120Diff_CellReg,meanbout2_2sec_30Same_CellReg,meanbout2_2sec_30Diff_CellReg,meanbout2_2sec_120Same_CellReg,meanbout2_2sec_120Diff_CellReg,meanboutHab1_2sec_30Same_CellReg,meanboutHab1_2sec_30Diff_CellReg,meanboutHab1_2sec_120Same_CellReg,meanboutHab1_2sec_120Diff_CellReg,meanboutHab2_2sec_30Same_CellReg,meanboutHab2_2sec_30Diff_CellReg,meanboutHab2_2sec_120Same_CellReg,meanboutHab2_2sec_120Diff_CellReg)));
maxboutmean2sec_CellReg = max(max(cat(1,meanbout1_2sec_30Same_CellReg,meanbout1_2sec_30Diff_CellReg,meanbout1_2sec_120Same_CellReg,meanbout1_2sec_120Diff_CellReg,meanbout2_2sec_30Same_CellReg,meanbout2_2sec_30Diff_CellReg,meanbout2_2sec_120Same_CellReg,meanbout2_2sec_120Diff_CellReg,meanboutHab1_2sec_30Same_CellReg,meanboutHab1_2sec_30Diff_CellReg,meanboutHab1_2sec_120Same_CellReg,meanboutHab1_2sec_120Diff_CellReg,meanboutHab2_2sec_30Same_CellReg,meanboutHab2_2sec_30Diff_CellReg,meanboutHab2_2sec_120Same_CellReg,meanboutHab2_2sec_120Diff_CellReg)));

%30 mins Diff 2 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_2sec_30Diff_CellReg(sortbout1_2sec_30Diff_CellReg,:))
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
colorbar
title('Trial1 2-seconds')
subplot(2,1,2)
plot(meanbout1_2sec_30Diff_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('30mins Novel Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_30Diff_2sec_CellReg')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_2sec_30Diff_CellReg(sortbout2_2sec_30Diff_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Trial2 2-seconds')
subplot(2,1,2)
plot(meanbout2_2sec_30Diff_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('30mins Novel Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_30Diff_2sec_CellReg')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_30Diff_CellReg(sortboutHab1_2sec_30Diff_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds')
subplot(2,1,2)
plot(meanboutHab1_2sec_30Diff_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('30mins Novel Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Diff_2sec_CellReg')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_30Diff_CellReg(sortboutHab2_2sec_30Diff_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds')
subplot(2,1,2)
plot(meanboutHab2_2sec_30Diff_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('30mins Novel Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Diff_2sec_CellReg')

%30 mins Same 2 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_2sec_30Same_CellReg(sortbout1_2sec_30Same_CellReg,:))
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
colorbar
title('Trial1 2-seconds')
subplot(2,1,2)
plot(meanbout1_2sec_30Same_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('30mins Familiar Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_30Same_2sec_CellReg')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_2sec_30Same_CellReg(sortbout2_2sec_30Same_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Trial2 2-seconds ')
subplot(2,1,2)
plot(meanbout2_2sec_30Same_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('30mins Familiar Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_30Same_2sec_CellReg')
% saveas(gcf,'CellFiringRateAcrossDaysAllMice120mins_Deconvolved.jpeg')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_30Same_CellReg(sortboutHab1_2sec_30Same_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds')
subplot(2,1,2)
plot(meanboutHab1_2sec_30Same_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('30mins Familiar Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Same_2sec_CellReg')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_30Same_CellReg(sortboutHab2_2sec_30Same_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds')
subplot(2,1,2)
plot(meanboutHab2_2sec_30Same_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('30mins Familiar Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Same_2sec_CellReg')


%120 mins Same 2 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_2sec_120Same_CellReg(sortbout1_2sec_120Same_CellReg,:))
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
colorbar
title('Trial1 2-seconds')
subplot(2,1,2)
plot(meanbout1_2sec_120Same_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('120mins Familiar Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_120Same_2sec_CellReg')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_2sec_120Same_CellReg(sortbout2_2sec_120Same_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Trial2 2-seconds')
subplot(2,1,2)
plot(meanbout2_2sec_120Same_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('120mins Familiar Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_120Same_2sec_CellReg')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_120Same_CellReg(sortboutHab1_2sec_120Same_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds')
subplot(2,1,2)
plot(meanboutHab1_2sec_120Same_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('120mins Familiar Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Same_2sec_CellReg')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_120Same_CellReg(sortboutHab2_2sec_120Same_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds')
subplot(2,1,2)
plot(meanboutHab2_2sec_120Same_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('120mins Familiar Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Same_2sec_CellReg')


%120 mins Diff 2 seconds
%trial1
figure
subplot(2,1,1)
imagesc(bout1_2sec_120Diff_CellReg(sortbout1_2sec_120Diff_CellReg,:))
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
colorbar
title('Trial1 2-seconds')
subplot(2,1,2)
plot(meanbout1_2sec_120Diff_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('120mins Novel Registerd Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial1_120Diff_2sec_CellReg')

%trial2
figure
subplot(2,1,1)
imagesc(bout2_2sec_120Diff_CellReg(sortbout2_2sec_120Diff_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Trial2 2-seconds')
subplot(2,1,2)
plot(meanbout2_2sec_120Diff_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('120mins Novel Registerd Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Trial2_120Diff_2sec_CellReg')

%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_120Diff_CellReg(sortboutHab1_2sec_120Diff_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds')
subplot(2,1,2)
plot(meanboutHab1_2sec_120Diff_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('120mins Novel Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Diff_2sec_CellReg')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_120Diff_CellReg(sortboutHab2_2sec_120Diff_CellReg,:))
caxis([minbout2sec_CellReg maxbout2sec_CellReg])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds')
subplot(2,1,2)
plot(meanboutHab2_2sec_120Diff_CellReg)
ylim([minboutmean2sec_CellReg maxboutmean2sec_CellReg])
sgtitle('120mins Novel Registered Cells')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Diff_2sec_CellReg')


%% Hab splithalf figures 

% First half 30 mins same
%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_30Same_1h(sortboutHab1_2sec_30Same_1h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds First Half')
subplot(2,1,2)
plot(meanboutHab1_2sec_30Same_1h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Same_2sec_1h')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_30Same_1h(sortboutHab2_2sec_30Same_1h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds First Half')
subplot(2,1,2)
plot(meanboutHab2_2sec_30Same_1h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Same_2sec_1h')

% Second half 30 mins same
%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_30Same_2h(sortboutHab1_2sec_30Same_2h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds Second Half')
subplot(2,1,2)
plot(meanboutHab1_2sec_30Same_2h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Same_2sec_2h')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_30Same_2h(sortboutHab2_2sec_30Same_2h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds Second Half')
subplot(2,1,2)
plot(meanboutHab2_2sec_30Same_2h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Same_2sec_2h')

%%%%%%


% First half 120 mins same
%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_120Same_1h(sortboutHab1_2sec_120Same_1h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds First Half')
subplot(2,1,2)
plot(meanboutHab1_2sec_120Same_1h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Same_2sec_1h')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_120Same_1h(sortboutHab2_2sec_120Same_1h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds First Half')
subplot(2,1,2)
plot(meanboutHab2_2sec_120Same_1h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Same_2sec_1h')

% Second half 120 mins same
%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_120Same_2h(sortboutHab1_2sec_120Same_2h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds Second Half')
subplot(2,1,2)
plot(meanboutHab1_2sec_120Same_2h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Same_2sec_2h')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_120Same_2h(sortboutHab2_2sec_120Same_2h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds Second Half')
subplot(2,1,2)
plot(meanboutHab2_2sec_120Same_2h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Familiar')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Same_2sec_2h')

%%%%%%%%%%%


% First half 30 mins Diff
%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_30Diff_1h(sortboutHab1_2sec_30Diff_1h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds First Half')
subplot(2,1,2)
plot(meanboutHab1_2sec_30Diff_1h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Diff_2sec_1h')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_30Diff_1h(sortboutHab2_2sec_30Diff_1h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds First Half')
subplot(2,1,2)
plot(meanboutHab2_2sec_30Diff_1h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Diff_2sec_1h')

% Second half 30 mins Diff
%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_30Diff_2h(sortboutHab1_2sec_30Diff_2h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds Second Half')
subplot(2,1,2)
plot(meanboutHab1_2sec_30Diff_2h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_30Diff_2sec_2h')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_30Diff_2h(sortboutHab2_2sec_30Diff_2h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds Second Half')
subplot(2,1,2)
plot(meanboutHab2_2sec_30Diff_2h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('30mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_30Diff_2sec_2h')

%%%%%%


% First half 120 mins Diff
%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_120Diff_1h(sortboutHab1_2sec_120Diff_1h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds First Half')
subplot(2,1,2)
plot(meanboutHab1_2sec_120Diff_1h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Diff_2sec_1h')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_120Diff_1h(sortboutHab2_2sec_120Diff_1h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds First Half')
subplot(2,1,2)
plot(meanboutHab2_2sec_120Diff_1h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Diff_2sec_1h')

% Second half 120 mins same
%Hab1
figure
subplot(2,1,1)
imagesc(boutHab1_2sec_120Diff_2h(sortboutHab1_2sec_120Diff_2h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 1 2-seconds Second Half')
subplot(2,1,2)
plot(meanboutHab1_2sec_120Diff_2h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab1_120Diff_2sec_2h')

%Hab2
figure
subplot(2,1,1)
imagesc(boutHab2_2sec_120Diff_2h(sortboutHab2_2sec_120Diff_2h,:))
caxis([minbout2sec maxbout2sec])
xticks([1:30:secframes2*2+1])
xticklabels({'-2', '-1', '0', '1', '2'})
colorbar
title('Habituation 2 2-seconds Second Half')
subplot(2,1,2)
plot(meanboutHab2_2sec_120Diff_2h)
ylim([minboutmean2sec maxboutmean2sec])
sgtitle('120mins Novel')
vline(secframes2,'k')
axis off
saveas(gcf,'Bouts_Hab2_120Diff_2sec_2h')

%% Passed Interaction percentage

x = [0,1];
y = x;

passed30minsSame1 = [];
passed30minsSame2 = [];
passed30minsSameh1 = [];
passed30minsSameh2 = [];
passed30minsDiff1 = [];
passed30minsDiff2 = [];
passed30minsDiffh1 = [];
passed30minsDiffh2 = [];
passed120minsSame1 = [];
passed120minsSame2 = [];
passed120minsSameh1 = [];
passed120minsSameh2 = [];
passed120minsDiff1 = [];
passed120minsDiff2 = [];
passed120minsDiffh1 = [];
passed120minsDiffh2 = [];

passed30minsSame1_2 = [];
passed30minsSame2_2 = [];
passed30minsSameh1_2 = [];
passed30minsSameh2_2 = [];
passed30minsDiff1_2 = [];
passed30minsDiff2_2 = [];
passed30minsDiffh1_2 = [];
passed30minsDiffh2_2 = [];
passed120minsSame1_2 = [];
passed120minsSame2_2 = [];
passed120minsSameh1_2 = [];
passed120minsSameh2_2 = [];
passed120minsDiff1_2 = [];
passed120minsDiff2_2 = [];
passed120minsDiffh1_2 = [];
passed120minsDiffh2_2 = [];

for i = 1 : length(perInteractionsZ(:,1,1))
    if ~isempty(perInteractionsZ{i,1,1})
        passed30minsSame1 = cat(2,passed30minsSame1,perInteractionsZ{i,1,1});
        passed30minsSame2 = cat(2,passed30minsSame2,perInteractionsZ{i,2,1});
        passed30minsSameh1 = cat(2,passed30minsSameh1,perInteractionsZ{i,3,1});
        passed30minsSameh2 = cat(2,passed30minsSameh2,perInteractionsZ{i,4,1});        
        
        passed30minsDiff1 = cat(2,passed30minsDiff1,perInteractionsZ{i,1,2});
        passed30minsDiff2 = cat(2,passed30minsDiff2,perInteractionsZ{i,2,2});
        passed30minsDiffh1 = cat(2,passed30minsDiffh1,perInteractionsZ{i,3,2});
        passed30minsDiffh2 = cat(2,passed30minsDiffh2,perInteractionsZ{i,4,2});
        
        passed30minsSame1_2 = cat(2,passed30minsSame1_2,perInteractionsZ_2{i,1,1});
        passed30minsSame2_2 = cat(2,passed30minsSame2_2,perInteractionsZ_2{i,2,1});
        passed30minsSameh1_2 = cat(2,passed30minsSameh1_2,perInteractionsZ_2{i,3,1});
        passed30minsSameh2_2 = cat(2,passed30minsSameh2_2,perInteractionsZ_2{i,4,1});        
        
        passed30minsDiff1_2 = cat(2,passed30minsDiff1_2,perInteractionsZ_2{i,1,2});
        passed30minsDiff2_2 = cat(2,passed30minsDiff2_2,perInteractionsZ_2{i,2,2});
        passed30minsDiffh1_2 = cat(2,passed30minsDiffh1_2,perInteractionsZ_2{i,3,2});
        passed30minsDiffh2_2 = cat(2,passed30minsDiffh2_2,perInteractionsZ_2{i,4,2});
    else
        passed120minsSame1 = cat(2,passed120minsSame1,perInteractionsZ{i,1,3});
        passed120minsSame2 = cat(2,passed120minsSame2,perInteractionsZ{i,2,3});
        passed120minsSameh1 = cat(2,passed120minsSameh1,perInteractionsZ{i,3,3});
        passed120minsSameh2 = cat(2,passed120minsSameh2,perInteractionsZ{i,4,3});
        
        passed120minsDiff1 = cat(2,passed120minsDiff1,perInteractionsZ{i,1,4});
        passed120minsDiff2 = cat(2,passed120minsDiff2,perInteractionsZ{i,2,4});
        passed120minsDiffh1 = cat(2,passed120minsDiffh1,perInteractionsZ{i,3,4});
        passed120minsDiffh2 = cat(2,passed120minsDiffh2,perInteractionsZ{i,4,4});
        
        passed120minsSame1_2 = cat(2,passed120minsSame1,perInteractionsZ{i,1,3});
        passed120minsSame2_2 = cat(2,passed120minsSame2,perInteractionsZ{i,2,3});
        passed120minsSameh1_2 = cat(2,passed120minsSameh1,perInteractionsZ{i,3,3});
        passed120minsSameh2_2 = cat(2,passed120minsSameh2,perInteractionsZ{i,4,3});
        
        passed120minsDiff1_2 = cat(2,passed120minsDiff1_2,perInteractionsZ_2{i,1,4});
        passed120minsDiff2_2 = cat(2,passed120minsDiff2_2,perInteractionsZ_2{i,2,4});
        passed120minsDiffh1_2 = cat(2,passed120minsDiffh1_2,perInteractionsZ_2{i,3,4});
        passed120minsDiffh2_2 = cat(2,passed120minsDiffh2_2,perInteractionsZ_2{i,4,4});
    end
end
m = 2;
n = 3;
c = 1;
%30 mins Same
subplot(m,n,c)
scatter(passed30minsSame1,passed30minsSame2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Trial 2 Familiar')
xlabel('Trial1')
ylabel('Trial2')
p = polyfit(passed30minsSame1,passed30minsSame2, 1);
px = [min(passed30minsSame1) max(passed30minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsSame1,passed30minsSameh1)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Hab1 Familiar')
xlabel('Trial1')
ylabel('Habituation1')
p = polyfit(passed30minsSame1,passed30minsSameh1, 1);
px = [min(passed30minsSame1) max(passed30minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsSame1,passed30minsSameh2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Hab2 Familiar')
xlabel('Trial1')
ylabel('Habituation2')
p = polyfit(passed30minsSame1,passed30minsSameh2, 1);
px = [min(passed30minsSame1) max(passed30minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsSameh1,passed30minsSame2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab1 vs Trial 2 Familiar')
xlabel('Habituation1')
ylabel('Trial2')
p = polyfit(passed30minsSameh1,passed30minsSame2, 1);
px = [min(passed30minsSame1) max(passed30minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsSameh2,passed30minsSame2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab2 vs Trial 2 Familiar')
xlabel('Habituation2')
ylabel('Trial2')
p = polyfit(passed30minsSameh2,passed30minsSame2, 1);
px = [min(passed30minsSame1) max(passed30minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsSameh1,passed30minsSameh2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab1 vs Hab1 Familiar')
xlabel('Habituation1')
ylabel('Habituation2')
p = polyfit(passed30minsSameh1,passed30minsSameh2, 1);
px = [min(passed30minsSame1) max(passed30minsSame1)];
py = polyval(p, px);
plot(px,py,'k')
sgtitle('Proportion of Interactions that passed Shuffling: 30 mins Familiar')
saveas(gcf,'PassedShuffledInteractionProportion_30minsSame')


%30 mins Novel
c = 1;
figure
subplot(m,n,c)
scatter(passed30minsDiff1,passed30minsDiff2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Trial 2 Novel')
xlabel('Trial1')
ylabel('Trial2')
p = polyfit(passed30minsDiff1,passed30minsDiff2, 1);
px = [min(passed30minsDiff1) max(passed30minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsDiff1,passed30minsDiffh1)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Hab1 Novel')
xlabel('Trial1')
ylabel('Habituation1')
p = polyfit(passed30minsDiff1,passed30minsDiffh1, 1);
px = [min(passed30minsDiff1) max(passed30minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsDiff1,passed30minsDiffh2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Hab2 Novel')
xlabel('Trial1')
ylabel('Habituation2')
p = polyfit(passed30minsDiff1,passed30minsDiff2, 1);
px = [min(passed30minsDiff1) max(passed30minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsDiffh1,passed30minsDiff2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab1 vs Trial 2 Novel')
xlabel('Habituation1')
ylabel('Trial2')
p = polyfit(passed30minsDiffh1,passed30minsDiff2, 1);
px = [min(passed30minsDiff1) max(passed30minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsDiffh2,passed30minsDiff2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab2 vs Trial 2 Novel')
xlabel('Habituation2')
ylabel('Trial2')
p = polyfit(passed30minsDiffh2,passed30minsDiff2, 1);
px = [min(passed30minsDiff1) max(passed30minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed30minsDiffh1,passed30minsDiffh2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab1 vs Hab1 Novel')
xlabel('Habituation1')
ylabel('Habituation2')
p = polyfit(passed30minsDiffh1,passed30minsDiffh2, 1);
px = [min(passed30minsDiff1) max(passed30minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')
sgtitle('Proportion of Interactions that passed Shuffling: 30 mins Novel')
saveas(gcf,'PassedShuffledInteractionProportion_30minsDiff')

%120 mins Same
c = 1;
figure
subplot(m,n,c)
scatter(passed120minsSame1,passed120minsSame2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Trial 2 Familiar')
xlabel('Trial1')
ylabel('Trial2')
p = polyfit(passed120minsSame1,passed120minsSame2, 1);
px = [min(passed120minsSame1) max(passed120minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsSame1,passed120minsSameh1)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Hab1 Familiar')
xlabel('Trial1')
ylabel('Habituation1')
p = polyfit(passed120minsSame1,passed120minsSameh1, 1);
px = [min(passed120minsSame1) max(passed120minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsSame1,passed120minsSameh2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Hab2 Familiar')
xlabel('Trial1')
ylabel('Habituation2')
p = polyfit(passed120minsSame1,passed120minsSameh2, 1);
px = [min(passed120minsSame1) max(passed120minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsSameh1,passed120minsSame2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab1 vs Trial 2 Familiar')
xlabel('Habituation1')
ylabel('Trial2')
p = polyfit(passed120minsSameh1,passed120minsSame2, 1);
px = [min(passed120minsSame1) max(passed120minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsSameh2,passed120minsSame2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab2 vs Trial 2 Familiar')
xlabel('Habituation2')
ylabel('Trial2')
p = polyfit(passed120minsSameh2,passed120minsSame2, 1);
px = [min(passed120minsSame1) max(passed120minsSame1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsSameh1,passed120minsSameh2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab1 vs Hab1 Familiar')
xlabel('Habituation1')
ylabel('Habituation2')
p = polyfit(passed120minsSameh1,passed120minsSameh2, 1);
px = [min(passed120minsSame1) max(passed120minsSame1)];
py = polyval(p, px);
plot(px,py,'k')
sgtitle('Proportion of Interactions that passed Shuffling: 120 mins Familiar')
saveas(gcf,'PassedShuffledInteractionProportion_120minsSame')


%120 mins Novel
c = 1;
figure
subplot(m,n,c)
scatter(passed120minsDiff1,passed120minsDiff2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Trial 2 Novel')
xlabel('Trial1')
ylabel('Trial2')
p = polyfit(passed120minsDiff1,passed120minsDiff2, 1);
px = [min(passed120minsDiff1) max(passed120minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsDiff1,passed120minsDiffh1)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Hab1 Novel')
xlabel('Trial1')
ylabel('Habituation1')
p = polyfit(passed120minsDiff1,passed120minsDiffh1, 1);
px = [min(passed120minsDiff1) max(passed120minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsDiff1,passed120minsDiffh2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Trial 1 vs Hab2 Novel')
xlabel('Trial1')
ylabel('Habituation2')
p = polyfit(passed120minsDiff1,passed120minsDiff2, 1);
px = [min(passed120minsDiff1) max(passed120minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsDiffh1,passed120minsDiff2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab1 vs Trial 2 Novel')
xlabel('Habituation1')
ylabel('Trial2')
p = polyfit(passed120minsDiffh1,passed120minsDiff2, 1);
px = [min(passed120minsDiff1) max(passed120minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsDiffh2,passed120minsDiff2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab2 vs Trial 2 Novel')
xlabel('Habituation2')
ylabel('Trial2')
p = polyfit(passed120minsDiffh2,passed120minsDiff2, 1);
px = [min(passed120minsDiff1) max(passed120minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')

c = c +1;
subplot(m,n,c)
scatter(passed120minsDiffh1,passed120minsDiffh2)
hold on
ylim([0 1])
xlim([0 1])
plot(x,y)
title('Hab1 vs Hab1 Novel')
xlabel('Habituation1')
ylabel('Habituation2')
p = polyfit(passed120minsDiffh1,passed120minsDiffh2, 1);
px = [min(passed120minsDiff1) max(passed120minsDiff1)];
py = polyval(p, px);
plot(px,py,'k')
sgtitle('Proportion of Interactions that passed Shuffling: 120 mins Novel')
saveas(gcf,'PassedShuffledInteractionProportion_120minsDiff')

%%
%{
passedcount1same_30 = 0;
passedcount2same_30 = 0;
passedcount1same_120 = 0;
passedcount2same_120 = 0;

passedcount1diff_30 = 0;
passedcount2diff_30 = 0;
passedcount1diff_120 = 0;
passedcount2diff_120 = 0;

passSame30 = 0;
passSame120 = 0;
passDiff30 = 0;
passDiff120 = 0;

cellcount = 0;
figure
hold on
for i = 1 : length(preInteractions(:,1,1))
    if ~isempty(preInteractions{i,1,1}) && ~isempty(preInteractions{i,2,2})
        scatter(preInteractions{i,1,1},preInteractions{i,2,1},'b')                 %Same 30 mins
        scatter(preInteractions{i,1,2},preInteractions{i,2,2},'g')                 %Different 30 mins
        
        passedcount1same_30 = passedcount1same_30 + length(find(preInteractions{i,1,1} > 0.5));
        passedcount1diff_30 = passedcount1diff_30 + length(find(preInteractions{i,1,2} > 0.5));
        passedcount2same_30 = passedcount2same_30 + length(find(preInteractions{i,2,1} > 0.5));
        passedcount2diff_30 = passedcount2diff_30 + length(find(preInteractions{i,2,2} > 0.5));
        
        passSame30 = passSame30 + length(intersect(find(preInteractions{i,1,1} > 0.5),find(preInteractions{i,2,1} > 0.5)));
        passDiff30 = passDiff30 + length(intersect(find(preInteractions{i,1,2} > 0.5),find(preInteractions{i,2,2} > 0.5)));
        
        cellcount = cellcount + length(preInteractions{i,1,1});
        cellcount = cellcount + length(preInteractions{i,1,2});
    elseif ~isempty(preInteractions{i,1,3}) && ~isempty(preInteractions{i,2,4})
        scatter(preInteractions{i,1,3},preInteractions{i,2,3},'k')                 %same 120mins
        scatter(preInteractions{i,1,4},preInteractions{i,2,4},'m')                 %different 120 mins
        
        passedcount1same_120 = passedcount1same_120 + length(find(preInteractions{i,1,3} > 0.5));
        passedcount1diff_120 = passedcount1diff_120 + length(find(preInteractions{i,1,4} > 0.5));
        passedcount2same_120 = passedcount2same_120 + length(find(preInteractions{i,2,3} > 0.5));
        passedcount2diff_120 = passedcount2diff_120 + length(find(preInteractions{i,2,4} > 0.5));
        
        passSame120 = passSame120 + length(intersect(find(preInteractions{i,1,3} > 0.5),find(preInteractions{i,2,3} > 0.5)));
        passDiff120 = passDiff120 + length(intersect(find(preInteractions{i,1,4} > 0.5),find(preInteractions{i,2,4} > 0.5)));
        
        cellcount = cellcount + length(preInteractions{i,1,3});
        cellcount = cellcount + length(preInteractions{i,1,4});
    end
end
legend({'Same 120 mins', 'Different 120 mins','Same 30 mins', 'Different 30 mins'})
ylabel('Proportion of Passed Interactions Trial1')
xlabel('Proportion of Passed Interactions Trial2')
title('Interactions Reliability Across all Cells')

figure
scatter([ones(length(find(sameZonemean1_30(:,1))))], [sameZonemean1_30(find(sameZonemean1_30(:,1)),1)],'b')
hold on
scatter([2*ones(length(find(sameZonemean2_30(:,1))))], [sameZonemean2_30(find(sameZonemean2_30(:,1)),1)],'b')
scatter([3*ones(length(find(diffZonemean1_30(:,1))))], [diffZonemean1_30(find(diffZonemean1_30(:,1)),1)],'g')
scatter([4*ones(length(find(diffZonemean2_30(:,1))))], [diffZonemean2_30(find(diffZonemean2_30(:,1)),1)],'g')

scatter([ones(length(find(sameZonemean1_120(:,1))))], [sameZonemean1_120(find(sameZonemean1_120(:,1)),1)],'m')
scatter([2*ones(length(find(sameZonemean2_120(:,1))))], [sameZonemean2_120(find(sameZonemean2_120(:,1)),1)],'m')
scatter([3*ones(length(find(diffZonemean1_120(:,1))))], [diffZonemean1_120(find(diffZonemean1_120(:,1)),1)],'r')
scatter([4*ones(length(find(diffZonemean2_120(:,1))))], [diffZonemean2_120(find(diffZonemean2_120(:,1)),1)],'r')
xlim([0.6 4.4])
xticks([1,2,3,4])
xticklabels({'Trial1 Same', 'Trial2 Same','Trial1 Different','Trial1 Different'})
title('Mean interactions Activity across animals')
%}


%% Functions

function [boutout] = boutcat(enterzones,exitzones, secframes,gcellind,zdata,boutout,i)
if i > 1 && enterzones(i) - exitzones(i-1) >= secframes
    boutout = cat(1,boutout,zdata(gcellind,enterzones(i)-secframes:enterzones(i)+secframes));
elseif i == 1
    boutout = cat(1,boutout,zdata(gcellind,enterzones(i)-secframes:enterzones(i)+secframes));
end
end

function [boutout] = boutcat_1h(enterzones,exitzones, secframes,gcellind,zdata,boutout,i)
if i > 1 && enterzones(i) - exitzones(i-1) >= secframes && enterzones(i) < round(length(zdata(1,:))/2)
    boutout = cat(1,boutout,zdata(gcellind,enterzones(i)-secframes:enterzones(i)+secframes));
elseif i == 1
    boutout = cat(1,boutout,zdata(gcellind,enterzones(i)-secframes:enterzones(i)+secframes));
end
end

function [boutout] = boutcat_2h(enterzones,exitzones, secframes,gcellind,zdata,boutout,i)
if i > 1 && enterzones(i) - exitzones(i-1) >= secframes && enterzones(i) > round(length(zdata(1,:))/2)
    boutout = cat(1,boutout,zdata(gcellind,enterzones(i)-secframes:enterzones(i)+secframes));
end
end

function [pass] = passCellReg(regind,ind_1reg,ind_2reg,ind_h1reg,ind_h2reg,ind_1reg_not,ind_2reg_not,ind_h1reg_not,ind_h2reg_not,ind_1regprev,ind_2regprev,ind_h1regprev,ind_h2regprev,ind_1regprev_not,ind_2regprev_not,ind_h1regprev_not,ind_h2regprev_not)
pass.pass_S1_D2 = regind(intersect(ind_1reg,ind_2regprev),:);
pass.pass_S1_Dn2 = regind(intersect(ind_1reg,ind_2regprev_not),:);
pass.pass_Sn1_D2= regind(intersect(ind_1reg_not,ind_2regprev),:);
pass.pass_Sn1_Dn2 = regind(intersect(ind_1reg_not,ind_2regprev_not),:);
pass.pass_S1_Dh2 = regind(intersect(ind_1reg,ind_h2regprev),:);
pass.pass_Sn1_Dh2 = regind(intersect(ind_1reg_not,ind_h2regprev),:);
pass.pass_S1_Dnh2 = regind(intersect(ind_1reg,ind_h2regprev_not),:);
pass.pass_Sn1_Dnh2 = regind(intersect(ind_1reg_not,ind_h2regprev_not),:);
pass.pass_Sh1_D2 = regind(intersect(ind_h1reg,ind_h2regprev),:);
pass.pass_Snh1_D2 = regind(intersect(ind_h1reg_not,ind_h2regprev),:);
pass.pass_Sh1_Dn2 = regind(intersect(ind_h1reg,ind_h2regprev_not),:);
pass.pass_Snh1_Dn2 = regind(intersect(ind_h1reg_not,ind_h2regprev_not),:);
pass.pass_S1_D1 = regind(intersect(ind_1reg,ind_1regprev),:);
pass.pass_S1_Dn1 = regind(intersect(ind_1reg,ind_1regprev_not),:);
pass.pass_Sn1_D1= regind(intersect(ind_1reg_not,ind_1regprev),:);
pass.pass_Sn1_Dn1 = regind(intersect(ind_1reg_not,ind_1regprev_not),:);
pass.pass_S2_D2 = regind(intersect(ind_2reg,ind_2regprev),:);
pass.pass_S2_Dn2 = regind(intersect(ind_2reg,ind_2regprev_not),:);
pass.pass_Sn2_D2= regind(intersect(ind_2reg_not,ind_2regprev),:);
pass.pass_Sn2_Dn2 = regind(intersect(ind_2reg_not,ind_2regprev_not),:);
pass.pass_S2_D1 = regind(intersect(ind_2reg,ind_1regprev),:);
pass.pass_S2_Dn1 = regind(intersect(ind_2reg,ind_1regprev_not),:);
pass.pass_Sn2_D1= regind(intersect(ind_2reg_not,ind_1regprev),:);
pass.pass_Sn2_Dn1 = regind(intersect(ind_2reg_not,ind_1regprev_not),:);
pass.pass_Sh2_Dh2 = regind(intersect(ind_h2reg,ind_h2regprev),:);
pass.pass_Sh2_Dnh2 = regind(intersect(ind_h2reg,ind_h2regprev_not),:);
pass.pass_Snh2_Dh2= regind(intersect(ind_h2reg_not,ind_h2regprev),:);
pass.pass_Snh2_Dnh2 = regind(intersect(ind_h2reg_not,ind_h2regprev_not),:);
pass.pass_Sh1_Dh1 = regind(intersect(ind_h1reg,ind_h1regprev),:);
pass.pass_Sh1_Dnh1 = regind(intersect(ind_h1reg,ind_h1regprev_not),:);
pass.pass_Snh1_Dh1= regind(intersect(ind_h1reg_not,ind_h1regprev),:);
pass.pass_Snh1_Dnh1 = regind(intersect(ind_h1reg_not,ind_h1regprev_not),:);
pass.pass_Sh2_Dh1 = regind(intersect(ind_h2reg,ind_h1regprev),:);
pass.pass_Sh2_Dnh1 = regind(intersect(ind_h2reg,ind_h1regprev_not),:);
pass.pass_Snh2_Dh1= regind(intersect(ind_h2reg_not,ind_h1regprev),:);
pass.pass_Snh2_Dnh1 = regind(intersect(ind_h2reg_not,ind_h1regprev_not),:);
pass.pass_Sh1_Dh2 = regind(intersect(ind_h1reg,ind_h2regprev),:);
pass.pass_Sh1_Dnh2 = regind(intersect(ind_h1reg,ind_h2regprev_not),:);
pass.pass_Snh1_Dh2= regind(intersect(ind_h1reg_not,ind_h2regprev),:);
pass.pass_Snh1_Dnh2 = regind(intersect(ind_h1reg_not,ind_h2regprev_not),:);
pass.pass_S1_S2 = regind(intersect(ind_1reg,ind_2reg),:);
pass.pass_S1_Sn2 = regind(intersect(ind_1reg,ind_2reg_not),:);
pass.pass_Sn1_S2= regind(intersect(ind_1reg_not,ind_2reg),:);
pass.pass_Sn1_Sn2 = regind(intersect(ind_1reg_not,ind_2reg_not),:);
pass.pass_S2_S2 = regind(intersect(ind_2reg,ind_2reg),:);
pass.pass_S2_Sn2 = regind(intersect(ind_2reg,ind_2reg_not),:);
pass.pass_Sn2_S2= regind(intersect(ind_2reg_not,ind_2reg),:);
pass.pass_Sn2_Sn2 = regind(intersect(ind_2reg_not,ind_2reg_not),:);
pass.pass_D1_D2 = regind(intersect(ind_1regprev,ind_2regprev),:);
pass.pass_D1_Dn2 = regind(intersect(ind_1regprev,ind_2regprev_not),:);
pass.pass_Dn1_D2= regind(intersect(ind_1regprev_not,ind_2regprev),:);
pass.pass_Dn1_Dn2 = regind(intersect(ind_1regprev_not,ind_2regprev_not),:);
pass.pass_D2_D2 = regind(intersect(ind_2regprev,ind_2regprev),:);
pass.pass_D2_Dn2 = regind(intersect(ind_2regprev,ind_2regprev_not),:);
pass.pass_Dn2_D2= regind(intersect(ind_2regprev_not,ind_2regprev),:);
pass.pass_Dn2_Dn2 = regind(intersect(ind_2regprev_not,ind_2regprev_not),:);
end

function [pass] = passNonCellReg_same(ind_1reg,ind_2reg,ind_h1reg,ind_h2reg,ind_1reg_not,ind_2reg_not,ind_h1reg_not,ind_h2reg_not)
pass.pass_S1_S2 = (intersect(ind_1reg,ind_2reg));
pass.pass_S1_Sn2 = (intersect(ind_1reg,ind_2reg_not));
pass.pass_Sn1_S2= (intersect(ind_1reg_not,ind_2reg));
pass.pass_Sn1_Sn2 = (intersect(ind_1reg_not,ind_2reg_not));
pass.pass_S2_S2 = (intersect(ind_2reg,ind_2reg));
pass.pass_S2_Sn2 = (intersect(ind_2reg,ind_2reg_not));
pass.pass_Sn2_S2= (intersect(ind_2reg_not,ind_2reg));
pass.pass_Sn2_Sn2 = (intersect(ind_2reg_not,ind_2reg_not));

pass.pass_Sh1_S2 = (intersect(ind_h1reg,ind_2reg));
pass.pass_Sh1_Sn2 = (intersect(ind_h1reg,ind_2reg_not));
pass.pass_Snh1_S2= (intersect(ind_h1reg_not,ind_2reg));
pass.pass_Snh1_Sn2 = (intersect(ind_h1reg_not,ind_2reg_not));
pass.pass_Sh2_S2 = (intersect(ind_h2reg,ind_2reg));
pass.pass_Sh2_Sn2 = (intersect(ind_h2reg,ind_2reg_not));
pass.pass_Snh2_S2= (intersect(ind_h2reg_not,ind_2reg));
pass.pass_Snh2_Sn2 = (intersect(ind_h2reg_not,ind_2reg_not));

pass.pass_S1_Sh2 = (intersect(ind_1reg,ind_h2reg));
pass.pass_S1_Snh2 = (intersect(ind_1reg,ind_h2reg_not));
pass.pass_Sn1_Sh2= (intersect(ind_1reg_not,ind_h2reg));
pass.pass_Sn1_Snh2 = (intersect(ind_1reg_not,ind_h2reg_not));
pass.pass_S2_Sh2 = (intersect(ind_2reg,ind_h2reg));
pass.pass_S2_Snh2 = (intersect(ind_2reg,ind_h2reg_not));
pass.pass_Sn2_Sh2= (intersect(ind_2reg_not,ind_h2reg));
pass.pass_Sn2_Snh2 = (intersect(ind_2reg_not,ind_h2reg_not));

pass.pass_Sh1_Sh2 = (intersect(ind_h1reg,ind_h2reg));
pass.pass_Sh1_Snh2 = (intersect(ind_h1reg,ind_h2reg_not));
pass.pass_Snh1_Sh2= (intersect(ind_h1reg_not,ind_h2reg));
pass.pass_Snh1_Snh2 = (intersect(ind_h1reg_not,ind_h2reg_not));
pass.pass_Sh2_Sh2 = (intersect(ind_h2reg,ind_h2reg));
pass.pass_Sh2_Snh2 = (intersect(ind_h2reg,ind_h2reg_not));
pass.pass_Snh2_Sh2= (intersect(ind_h2reg_not,ind_h2reg));
pass.pass_Shn2_Shn2 = (intersect(ind_h2reg_not,ind_h2reg_not));
end

function [pass] = passNonCellReg_diff(ind_1reg,ind_2reg,ind_h1reg,ind_h2reg,ind_1reg_not,ind_2reg_not,ind_h1reg_not,ind_h2reg_not)
pass.pass_D1_D2 = (intersect(ind_1reg,ind_2reg));
pass.pass_D1_Dn2 = (intersect(ind_1reg,ind_2reg_not));
pass.pass_Dn1_D2= (intersect(ind_1reg_not,ind_2reg));
pass.pass_Dn1_Dn2 = (intersect(ind_1reg_not,ind_2reg_not));
pass.pass_D2_D2 = (intersect(ind_2reg,ind_2reg));
pass.pass_D2_Dn2 = (intersect(ind_2reg,ind_2reg_not));
pass.pass_Dn2_D2= (intersect(ind_2reg_not,ind_2reg));
pass.pass_Dn2_Dn2 = (intersect(ind_2reg_not,ind_2reg_not));

pass.pass_Dh1_D2 = (intersect(ind_h1reg,ind_2reg));
pass.pass_Dh1_Dn2 = (intersect(ind_h1reg,ind_2reg_not));
pass.pass_Dnh1_D2= (intersect(ind_h1reg_not,ind_2reg));
pass.pass_Dnh1_Dn2 = (intersect(ind_h1reg_not,ind_2reg_not));
pass.pass_Dh2_D2 = (intersect(ind_h2reg,ind_2reg));
pass.pass_Dh2_Dn2 = (intersect(ind_h2reg,ind_2reg_not));
pass.pass_Dnh2_D2= (intersect(ind_h2reg_not,ind_2reg));
pass.pass_Dnh2_Dn2 = (intersect(ind_h2reg_not,ind_2reg_not));

pass.pass_D1_Dh2 = (intersect(ind_1reg,ind_h2reg));
pass.pass_D1_Dnh2 = (intersect(ind_1reg,ind_h2reg_not));
pass.pass_Dn1_Dh2= (intersect(ind_1reg_not,ind_h2reg));
pass.pass_Dn1_Dnh2 = (intersect(ind_1reg_not,ind_h2reg_not));
pass.pass_D2_Dh2 = (intersect(ind_2reg,ind_h2reg));
pass.pass_D2_Dnh2 = (intersect(ind_2reg,ind_h2reg_not));
pass.pass_Dn2_Dh2= (intersect(ind_2reg_not,ind_h2reg));
pass.pass_Dn2_Dnh2 = (intersect(ind_2reg_not,ind_h2reg_not));

pass.pass_Dh1_Dh2 = (intersect(ind_h1reg,ind_h2reg));
pass.pass_Dh1_Dnh2 = (intersect(ind_h1reg,ind_h2reg_not));
pass.pass_Dnh1_Dh2= (intersect(ind_h1reg_not,ind_h2reg));
pass.pass_Dnh1_Dnh2 = (intersect(ind_h1reg_not,ind_h2reg_not));
pass.pass_Dh2_Dh2 = (intersect(ind_h2reg,ind_h2reg));
pass.pass_Dh2_Dnh2 = (intersect(ind_h2reg,ind_h2reg_not));
pass.pass_Dnh2_Dh2= (intersect(ind_h2reg_not,ind_h2reg));
pass.pass_Dhn2_Dhn2 = (intersect(ind_h2reg_not,ind_h2reg_not));
end