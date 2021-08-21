%% Create figure for concactenated Sessions SIT data

SITseshnum = 28;

mousename = SITproximityData.micenames{SITseshnum};

folderparts = strsplit(SITproximityData.folderpaths{SITseshnum},'\');
dateparts = strsplit(folderparts{end-3},' ');
seshdate = dateparts{2};

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
trial2 = ms.deconvolvedSig(startframe2+l1:timeframes+startframe2+l1,:)';

% hab1 =  normalize(hab1,2,'range');
% trial1 =  normalize(trial1,2,'range');
% hab2 =  normalize(hab2,2,'range');
% trail2 =  normalize(trial2,2,'range');

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
    
    if i == 1
        mult = 0;
        vline(enterzonesHab,'g:')
        vline(exitzonesHab,'r:')
    else
        mult = 1;
        vline(timeframes*(i)+enterzonesHab,'g:')
        vline(timeframes*(i)+exitzones,'r:')
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
            if ~isempty(hab2)
                mapframe(round(i/gridfactor),round(j/gridfactor)) = length(hab2);
                mapframeHab2(round(i/gridfactor),round(j/gridfactor)) = length(hab2);
            else
                nanmaskHab2(round(i/gridfactor),round(j/gridfactor)) = nan;
            end
            
%             if ~isempty(t1) && isempty(t2)
%                 nanmask2(round(i/gridfactor),round(j/gridfactor)) = nan;
%                 mapframe(round(i/gridfactor),round(j/gridfactor)) = length(t1);
%                 mapframe1(round(i/gridfactor),round(j/gridfactor)) = length(t1);
%             elseif isempty(t1) && ~isempty(t2)
%                 mapframe(round(i/gridfactor),round(j/gridfactor)) = length(t2);
%                 nanmask1(round(i/gridfactor),round(j/gridfactor)) = nan;
%                 mapframe2(round(i/gridfactor),round(j/gridfactor)) = length(t2);
%             else
%                 mapframe(round(i/gridfactor),round(j/gridfactor)) = length(t1) + length(t2);
%                 mapframe1(round(i/gridfactor),round(j/gridfactor)) = length(t1);
%                 mapframe2(round(i/gridfactor),round(j/gridfactor)) = length(t2);
%             end
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

% [lf1x, lf1y] = find(mapframe1 <= 9);
% [lf2x, lf2y]= find(mapframe2 <= 9);
% map1(lf1x,lf1y) = 0;
% map2(lf2x,lf2y) = 0;


filtkernel = 1;
count = 1;
%{
track = cat(1,track1,track2);
figure('Position', get(0, 'Screensize')); 
for i = 1 : celltotal    
    map(:,:,i) = imgaussfilt(map(:,:,i),filtkernel).*nanmask;    
    subplot(4,4,count)
    imagesc(map(:,:,i),'AlphaData',~isnan(nanmask))
    title(['Cell ' num2str(gcellind(i))])
    subplot(4,4,count+4)    
    plot(track(:,1),track(:,2))
    hold on
    indfiring = find(traces(i,:));
    scatter(track(indfiring,1),track(indfiring,2),'r.')
    set(gca,'visible','off')
    if count == 4
        count = 9;
    else
        count = count+1;
    end
    if count > 12
        count = 1;
        saveas(gcf,['FiringRateMap_' num2str(gcellind(i))])
        saveas(gcf,['FiringRateMap_' num2str(gcellind(i)) '.jpeg'])  
        pause(0.01)
        clf
    end
end
%}
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

%% Mean Interaction Activity
% ztrial = traces
ztrials = normalize(cat(2,hab1,trial1,hab2,trial2),2);
zhab1 = ztrials(:,1:length(hab1(1,:)));
ztrial1 = ztrials(:,length(zhab1(1,:))+1:length(zhab1(1,:))+length(trial1(1,:)));
zhab2 = ztrials(:,length(traces1(1,:))+1:length(traces1(1,:))+length(hab2(1,:)));
ztrial2 = ztrials(:,length(traces1(1,:))+length(hab2(1,:))+1:end);

if length(enterzones1) ~= length(exitzones1)
    if length(enterzones1) > length(exitzones1)
        enterzones1(end) = [];
    else
        exitzones1(end) = [];
    end
end
if length(enterzones2) ~= length(exitzones2)
    if length(enterzones2) > length(exitzones2)
        enterzones2(end) = [];
    else
        exitzones2(end) = [];
    end
end
if length(enterzonesHab1) ~= length(exitzonesHab1)
    if length(enterzonesHab1) > length(exitzonesHab1)
        enterzonesHab1(end) = [];
    else
        exitzonesHab1(end) = [];
    end
end
if length(enterzonesHab2) ~= length(exitzonesHab2)
    if length(enterzonesHab2) > length(exitzonesHab2)
        enterzonesHab2(end) = [];
    else
        exitzonesHab2(end) = [];
    end
end

timebinnum = 5; 
interactionsMean = zeros(celltotal,timebinnum);
zoneActivitysplit1 = zeros(celltotal,length(enterzones1),timebinnum);
for time = 1 : timebinnum
    for i = 1 : length(enterzones1)        
        binsize = round((exitzones1(i)-enterzones1(i))/timebinnum);
        if binsize > 1
            for j = 1 : celltotal
                zoneActivitysplit1(j,i,time) = sum(ztrial1(gcellind(j),enterzones1(i)+binsize*(time-1):enterzones1(i)+binsize*(time)))/binsize;
                zoneActivityRaw1(j,i) = sum(trial1(gcellind(j),enterzones1(i):exitzones1(i)))/(exitzones1(i)-enterzones1(i));
            end
        end
    end    
end
zoneActivitysplitHab1 = zeros(celltotal,length(enterzonesHab1),timebinnum);
for time = 1 : timebinnum
    for i = 1 : length(enterzonesHab1)        
        binsize = round((exitzonesHab1(i)-enterzonesHab1(i))/timebinnum);
        if binsize > 1
            for j = 1 : celltotal
                zoneActivitysplitHab1(j,i,time) = sum(zhab1(gcellind(j),enterzonesHab1(i)+binsize*(time-1):enterzonesHab1(i)+binsize*(time)))/binsize;
                zoneActivityHabRaw1(j,i) = sum(hab1(gcellind(j),enterzonesHab1(i):exitzonesHab1(i)))/(exitzonesHab1(i)-enterzonesHab1(i));
            end
        end
    end    
end
zoneActivitysplit2 = zeros(celltotal,length(enterzones2),timebinnum);
for time = 1 : timebinnum
    for i = 1 : length(enterzones2)        
        binsize = round((exitzones2(i)-enterzones2(i))/timebinnum);
        if binsize > 1
            for j = 1 : celltotal
                zoneActivitysplit2(j,i,time) = sum(ztrial2(gcellind(j),enterzones2(i)+binsize*(time-1):enterzones2(i)+binsize*(time)))/binsize;
                zoneActivityRaw2(j,i) = sum(trial2(gcellind(j),enterzones2(i):exitzones2(i)))/(exitzones2(i)-enterzones2(i));
            end
        end
    end
end
zoneActivitysplitHab2 = zeros(celltotal,length(enterzonesHab2),timebinnum);
for time = 1 : timebinnum
    for i = 1 : length(enterzonesHab2)
        binsize = round((exitzonesHab2(i)-enterzonesHab2(i))/timebinnum);
        if binsize > 1
            for j = 1 : celltotal
                zoneActivitysplitHab2(j,i,time) = sum(zhab2(gcellind(j),enterzonesHab2(i)+binsize*(time-1):enterzonesHab2(i)+binsize*(time)))/binsize;
                zoneActivityHabRaw2(j,i) = sum(hab2(gcellind(j),enterzonesHab2(i):exitzonesHab2(i)))/(exitzonesHab2(i)-enterzonesHab2(i));
            end
        end
    end
end


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

saveas(gcf,['MeanInteractionFiringRate'])
saveas(gcf,['MeanInteractionFiringRate.jpeg'])
