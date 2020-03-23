% Compile Signal to Noise ratio of miniscope recordings
% Script designed for looking at single mouse, across sessions
% Change mouse number for each 'save' command

% Whole session
% % % p = getFilePaths(pwd,'.mat');
snrWhole = cell(size(p,1),1);
for i = 1:size(p,1)
    s = load(p{i},'processed');
    j = s.processed.snr.whole;
    snrWhole{i,1} = j; 
end
save('snrWhole_4DH5.mat', 'snrWhole'); % Adjust accordingly
m = median(cell2mat(snrWhole));
figure, histogram(cell2mat(snrWhole)), hold on
plot([m, m], ylim, 'Color', 'b', 'LineWidth', 2);
text(35,350,'Median', 'Color', 'k', 'FontSize', 20)
text(35,300,num2str(m), 'Color', 'b', 'FontSize', 16)

% Get session names
ss = size(p,1);
session = string(NaN(ss,1));
for i = 1:ss
f = load(p{i}, 'ms');
r = string(f.ms.dirName);
session(i,1) = r;
end

% Array of whole SNR medians for each mouse session
snrWholeMeds = NaN(ss,1);
for k = 1:ss
    a = median(snrWhole{k,1});
    snrWholeMeds(k,1) = a;
end
save('snrWholeMeds_4DH5.mat', 'snrWholeMeds'); % Adjust accordingly

% Figures of whole SNR for each mouse session
for i = 1:size(p,1)
figure, histogram(snrWhole{i,1},0:2:60), hold on % Adjust accordingly
b = median(snrWhole{i,1});
plot([b, b], ylim, 'Color', 'b', 'LineWidth', 2)
title(num2str(p{i}))
text(22,30,'Median', 'Color', 'k', 'FontSize', 14)
text(22,26,num2str(b), 'Color', 'b', 'FontSize', 12)
end

% Split session
snrSplit = cell(size(p,1),2); 
for i = 1:size(p,1)
    s = load(p{i},'processed');
    j = s.processed.snr.splithalf(:,1);
    snrSplit{i,1} = j;
    k = s.processed.snr.splithalf(:,2);
    snrSplit{i,2} = k;
end
save('snrSplit_4DH5.mat', 'snrSplit'); % Adjust accordingly
n = median(cell2mat(snrSplit));
g = cell2mat(snrSplit);
figure, histogram(g(:,1), 0:2:60), hold on
histogram(g(:,2), 0:2:60)
plot([n(1), n(1)], ylim, 'Color', 'b', 'LineWidth', 2)
plot([n(2), n(2)], ylim, 'Color', 'r', 'LineWidth', 2)
text(35,500,'Medians', 'Color', 'k', 'FontSize', 20)
text(35,425,num2str(n(1)), 'Color', 'b', 'FontSize', 16)
text(35,350,num2str(n(2)), 'Color', 'r', 'FontSize', 16)

% Array of whole SNR medians for each mouse session
ss = size(p,1);
snrSplitMeds = NaN([ss 2]);
for k = 1:ss
    a = median(snrSplit{k,1});
    snrSplitMeds(k,1) = a;
    a = median(snrSplit{k,2});
    snrSplitMeds(k,2) = a;
end
save('snrSplitMeds_4DH5.mat', 'snrSplitMeds'); % Adjust accordingly

% Figures of first and second half SNR for each mouse session
for i = 1:size(p,1)
figure, histogram(snrSplit{i,1},0:2:60), hold on % Adjust accordingly
histogram(snrSplit{i,2},0:2:60)
b = median(snrSplit{i,1});
c = median(snrSplit{i,2});
plot([b, b], ylim, 'Color', 'b', 'LineWidth', 2)
plot([c, c], ylim, 'Color', 'r', 'LineWidth', 2)
title(num2str(p{i}))
text(25,35,'Medians', 'Color', 'k', 'FontSize', 14)
text(25,30,num2str(b), 'Color', 'b', 'FontSize', 12)
text(25,25,num2str(c), 'Color', 'r', 'FontSize', 12)
end
