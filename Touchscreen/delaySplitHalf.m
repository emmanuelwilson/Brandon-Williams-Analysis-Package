%%Splithalf comparison of delay periods
function [d1, d2] = delaySplitHalf(ms,dctrace, ditrace,dcctrace,dicctrace,dc,di,dcc,dic, delay, delay2) 
%Find the halfway seperation point, making sure it doesn't split a delay period
half = round(length(ms.FiltTraces(:,1))/2);
mind = 0;

% for i = 1: length(delay)
%     if i >= length(delay2)
%         break
%     elseif half>= delay(i) && half<delay2(i)
%         [~,mind] = min(abs(delay-half));
%         half = delay(mind)-1;
%     end
% end

if mind == 0 
    temp = double((delay-half)>0);
    temp(temp == 0) = 100000;
    [~,mind] = min(temp);
    half = delay(mind)-1;
end

dctrace1 = nan(length(dctrace(:,1,1)),length(dctrace(1,:,1)),length(dc));
ditrace1 = nan(length(dctrace(:,1,1)),length(ditrace(1,:,1)),length(di));
dcctrace1 = nan(length(dctrace(:,1,1)),length(dcctrace(1,:,1)),length(dcc));
dicctrace1 = nan(length(dctrace(:,1,1)),length(dicctrace(1,:,1)),length(dic));

dctrace2 = nan(length(dctrace(:,1,1)),length(dctrace(1,:,1)),length(dc));
ditrace2 = nan(length(dctrace(:,1,1)),length(ditrace(1,:,1)),length(di));
dcctrace2 = nan(length(dctrace(:,1,1)),length(dcctrace(1,:,1)),length(dcc));
dicctrace2 = nan(length(dctrace(:,1,1)),length(dicctrace(1,:,1)),length(dic));

dccount1 = 1;
dicount1 = 1;
dcccount1 = 1;
dicccount1 = 1;

dccount2 = 1;
dicount2 = 1;
dcccount2 = 1;
dicccount2 = 1;
    
for i =1: length(delay)
    if i < mind
        if max(delay(i) == dc) == 1
            dctrace1(:,:,dccount1) = dctrace(:,:,find(dc == delay(i)));
            dccount1 = dccount1 + 1;
        elseif max(delay(i) == di)==1
            ditrace1(:,:,dicount1) = ditrace(:,:,find(di == delay(i)));
            dicount1 = dicount1 + 1;
        elseif max(delay(i) == dcc)==1
            dcctrace1(:,:,dcccount1) = dcctrace(:,:,find(dcc == delay(i)));
            dcccount1 = dcccount1 + 1;
        elseif max(delay(i)== dic)==1
            dicctrace1(:,:,dicccount1) = dicctrace(:,:,find(dic == delay(i)));
            dicccount1 = dicccount1 + 1;
        end
    else
        if max(delay(i) == dc)==1
            dctrace2(:,:,dccount2) = dctrace(:,:,find(dc == delay(i)));
            dccount2 = dccount2 + 1;
        elseif max(delay(i) == di)==1
            ditrace2(:,:,dicount2) = ditrace(:,:,find(di == delay(i)));
            dicount2 = dicount2 + 1;
        elseif max(delay(i) == dcc)==1
            dcctrace2(:,:,dcccount2) = dcctrace(:,:,find(dcc == delay(i)));
            dcccount2 = dcccount2 + 1;
        elseif max(delay(i) == dic)==1
            dicctrace2(:,:,dicccount2) = dicctrace(:,:,find(dic == delay(i)));
            dicccount2 = dicccount2 + 1;
        end
    end
end

for i = 1 : 4
    if i == 1
        trace = dctrace;
    elseif i == 2
        trace = ictrace;
    elseif i == 3
        trace = dcctrace;
    elseif i ==4 
        trace = dictrace;
    end
    precorr = permute(trace,[2 3 1]);
    Corr = corrcoef(precorr(:,:,1));
    figure
    c(i) = imagesc(Corr)
    vline(mind,'r')
    hline(mind,'r')
    xlabel('Trial Number')
    ylabel('Trial Number')
    if i == 1
        title('Correct Trial Delay Correlation')
    elseif i == 2
        title('Incorrect Trial Delay Correlation')
    elseif i == 3
        title('Correct Correction Trial Delay Correlation')
    elseif i ==4 
        title('Incorrect Correction Trial Delay Correlation')
    end
end
savefig('Delay Correlation')

dcrast1 = nanmean(dctrace1,3);
dirast1 = nanmean(ditrace1,3);
dccrast1 = nanmean(dcctrace1,3);
dicrast1 = nanmean(dicctrace1,3);

%Zero the minimum values Delay periods
for i = 1 : length(dctrace(:,1,1))
    dcrast1(i,:) = dcrast1(i,:)-min(dcrast1(i,:));
end
for i = 1 : length(ditrace(:,1,1))
    dirast1(i,:) = dirast1(i,:)-min(dirast1(i,:));
end
for i = 1 : length(dcctrace(:,1,1))
    dccrast1(i,:) = dccrast1(i,:)-min(dccrast1(i,:));
end
for i = 1 : length(dicctrace(:,1,1))
    dicrast1(i,:) = dicrast1(i,:)-min(dicrast1(i,:));
end

dcrast1 = sortpeaks(dcrast1);
dirast1 = sortpeaks(dirast1);
dccrast1 = sortpeaks(dccrast1);
dicrast1 = sortpeaks(dicrast1);

%display the delay period
d(1) = figure;
h = imagesc(dcrast1);
% set(h,'LineStyle','none')
axis('tight')
colormap('hot')
colorbar
view(2)
xlabel('Time(frame)')
ylabel('Neuron Number')
title('Correct Trails First Half Delay Period Calcium Activity')
pause(0.01)

d(2) = figure;
h = imagesc(dirast1);
% set(h,'LineStyle','none')
axis('tight')
colormap('hot')
colorbar
view(2)
xlabel('Time(frame)')
ylabel('Neuron Number')
title('Incorrect Trails First Half Delay Period Calcium Activity')
pause(0.01)

d(3) = figure;
h = imagesc(dccrast1);
% set(h,'LineStyle','none')
axis('tight')
colormap('hot')
colorbar
view(2)
xlabel('Time(frame)')
ylabel('Neuron Number')
title('Correct Correction Trails First Half Delay Period Calcium Activity')
pause(0.01)

d(4) = figure;
h = imagesc(dicrast1);
% set(h,'LineStyle','none')
axis('tight')
colormap('hot')
colorbar
view(2)
xlabel('Time(frame)')
ylabel('Neuron Number')
title('Incorrect Correction Trails First Half Delay Period Calcium Activity')
pause(0.01)

savefig(d,'DelayFirstHalf')

dcrast2 = nanmean(dctrace2,3);
dirast2 = nanmean(ditrace2,3);
dccrast2 = nanmean(dcctrace2,3);
dicrast2 = nanmean(dicctrace2,3);

%Zero the minimum values Delay periods
for i = 1 : length(dctrace(:,1,1))
    dcrast2(i,:) = dcrast2(i,:)-min(dcrast2(i,:));
end
for i = 1 : length(ditrace(:,1,1))
    dirast2(i,:) = dirast2(i,:)-min(dirast2(i,:));
end
for i = 1 : length(dcctrace(:,1,1))
    dccrast2(i,:) = dccrast2(i,:)-min(dccrast2(i,:));
end
for i = 1 : length(dicctrace(:,1,1))
    dicrast2(i,:) = dicrast2(i,:)-min(dicrast2(i,:));
end

dcrast2 = sortpeaks(dcrast2);
dirast2 = sortpeaks(dirast2);
dccrast2 = sortpeaks(dccrast2);
dicrast2 = sortpeaks(dicrast2);

%display the delay period
d(1) = figure;
h = imagesc(dcrast2);
% set(h,'LineStyle','none')
axis('tight')
colormap('hot')
colorbar
view(2)
xlabel('Time(frame)')
ylabel('Neuron Number')
title('Correct Trails Second Half Delay Period Calcium Activity')
pause(0.01)

d(2) = figure;
h = imagesc(dirast2);
% set(h,'LineStyle','none')
axis('tight')
colormap('hot')
colorbar
view(2)
xlabel('Time(frame)')
ylabel('Neuron Number')
title('Incorrect Trails Second Half Delay Period Calcium Activity')
pause(0.01)

d(3) = figure;
h = imagesc(dccrast2);
% set(h,'LineStyle','none')
axis('tight')
colormap('hot')
colorbar
view(2)
xlabel('Time(frame)')
ylabel('Neuron Number')
title('Correct Correction Trails Second Half Delay Period Calcium Activity')
pause(0.01)

d(4) = figure;
h = imagesc(dicrast2);
% set(h,'LineStyle','none')
axis('tight')
colormap('hot')
colorbar
view(2)
xlabel('Time(frame)')
ylabel('Neuron Number')
title('Incorrect Correction Trails Second Half Delay Period Calcium Activity')
pause(0.01)

savefig(d,'DelaySecondHalf')

d1.dctrace = dctrace1;
d1.ditrace = ditrace1;
d1.dcctrace = dcctrace1;
d1.dicctrace = dicctrace1;

d2.dctrace = dctrace2;
d2.ditrace = ditrace2;
d2.dcctrace = dcctrace2;
d2.dicctrace = dicctrace2;

end

function rast = sortpeaks(rast)
[~,maxind] = max(rast,[],2);
[~, rastsort] = sort(maxind);
% rastsort = flipud(rastsort);
rast = rast(rastsort,:);
end