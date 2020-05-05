function ms = msRun(doPath)

%% Generates the initial ms data struct for data set contained in current folder

% cleanCalciumVideos(doPath);

ms = msGenerateVideoObj(doPath,'msCam');  

% Add the following line if analyzing a subset of data
%ms.time = linspace(0,ms.numFrames/30*1000,ms.numFrames);
% ------------------------------------------------------


%ms = msColumnCorrection_1(ms,5); %Generally not used

plotting = false;
ms = msSelectROIs(ms); %prompts the use to chose a frame in a video. 
ms = msAlignmentFFT(ms,plotting); %takes chosen frames and aligns them. 
downsample = 5;
ms = msMeanFrame_1(ms,downsample);
ms = msSelectAlignment(ms);
ms = msSelectMask(ms,downsample);  

%ms = msFluorFrameProps(ms);

%% Select fluorescnece thesh for good frames
%ms = msSelectFluorThresh(ms);

%% Allows user to select ROIs for each data folder
%%%%%%%ms = msSelectROIs(ms);
%% Run alignment across all ROIs
%%%%%%%plotting = true;
%%%%%%%tic
%%%%%%%ms = msAlignmentFFT(ms,plotting);
%%%%%%%toc
%% Calculate mean frames
%%%%%%%downsample = 5;
%%%%%%%ms = msMeanFrame(ms,downsample);

%% Manually inspect and select best alignment
%%%%%%%ms = msSelectAlignment(ms);

% Save video - Zaki
% msSaveVideo(ms,[], 1, true, true, true, 'savedVid1');


%% Segment Sessions
plotting = true;
ms = msFindBrightSpots(ms,30,[],.03,0.02,plotting);
ms = msAutoSegment2(ms,[],[60 700],downsample,.85,plotting);  % [60 700],5,.90

%% Calculate Segment relationships
calcCorr = false;
calcDist = true;
calcOverlap = true;
ms = msCalcSegmentRelations(ms, calcCorr, calcDist, calcOverlap);

%% Clean Segments
corrThresh = [];
distThresh = 7;
overlapThresh = .8;
ms = msCleanSegments(ms,corrThresh,distThresh,overlapThresh);

%% Calculate Segment relationships
calcCorr = false;
calcDist = true;
calcOverlap = true;
ms = msCalcSegmentRelations(ms, calcCorr, calcDist, calcOverlap);

%% Calculate segment centroids
ms = msSegmentCentroids(ms);

%% Extract dF/F
ms = msExtractdFFTraces(ms);
ms = msCleandFFTraces(ms);
ms = msExtractFiring(ms);

%% Align across sessions
% ms = msAlignBetweenSessions(msRef,ms);

%% Count segments in common field
% msBatchSegmentsInField(pwd);

%% Match segments across sessions
% distThresh = 5;
% msBatchMatchSegmentsBetweenSessions(pwd, distThresh);


%% BEHAV STUFF
% 
% %% Generate behav.m
% behav = msGenerateVideoObj(pwd,'behavCam');
% 
% %% Select ROI and HSV for tracking
% behav = msSelectPropsForTracking(behav); 
% 
% %% Extract position
% trackLength = 200;%cm
% behav = msExtractBehavoir(behav, trackLength); 



% open('HD_180_deg.fig');
% h = gcf; %current figure handle
%       axesObjs = get(h, 'Children');  %axes handles
%       dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
%       objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% ydata = get(dataObjs, 'YData');
% figure;plot(ms.trace(:,1))
% hold on
% plot(xdata, (ydata/180))


% figure;plot(ms.trace(:,1)*8)
%  hold on
%  plot(HD_180_deg/180)

end