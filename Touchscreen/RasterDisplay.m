%Population Level activity Raster plot display
function [] = RasterDisplay(rasters,FigName)

% mkdir(folder_name);
totrast = nanmean(rasters,3);
for i = 1: length(totrast(:,1))
    totrast(i,:) = normalize(totrast(i,:),'range');
end

for i =length(totrast(:,1)):-1:1
    m = mean(totrast(i,:));
    if ~((max(totrast(i,:))-m)>= (2*std(totrast(i,:))))
        totrast(i,:) = [];
    end
end

[~,maxind] = max(totrast,[],2);
[~, rastsort] = sort(maxind);
totrast = totrast(rastsort,:);
figure
imagesc(totrast)
title(FigName)
ylabel('Neuron Number')
xlabel('Time(frames)')
colormap(parula)
colorbar
savefig([FigName])

end