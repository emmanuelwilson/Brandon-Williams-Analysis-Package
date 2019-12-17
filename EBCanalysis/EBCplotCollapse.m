dist = sum(out.rm(:,:,7),2);
angle = sum(out.rm(:,:,7));
dist = zscore(dist);
angle = zscore(angle);
[~, distMax] = find(max(dist));
distThresh =0.2*(max(dist)-min(dist));
angleMax = find(max(angle));
angleThresh = 0.2*(max(angle)- min(angle));
distB = find(dist>=distThresh);

%{
% distA = find(dist>=1.65);
% angleA = find(angle>=1.65);
% p = find(diff(angleA) == 1);
% q = [p , max(p) + 1];
% angleB = angleA(q);
% p = find(diff(distA) == 1);
% q = [p' max(p) + 1];
% distB = distA(q);
%}
figure
plot(angle)
vline(angleB(1),'r')
vline(angleB(end),'r')
figure
plot(dist)
view(-90,-90)
vline(distB(1),'r')
vline(distB(end),'r')
figure
imagesc(out.rm(:,:,7))
vline(angleB(1),'w')
vline(angleB(end),'w')
hline(distB(1),'w')
hline(distB(end),'w')
max(dist)
max(angle)