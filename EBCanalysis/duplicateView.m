%% Takes your miniscope structure and list of potential duplicates to display

function [] = duplicateView(ms,dup)
traces = ms.FiltTraces;
footprints = ms.SFPs;
indOut = 1:length(ms.FiltTraces(1,:));
indOut(dup) = [];

traces(:,indOut) = [];
footprints(:,:,indOut) = [];

figure
for i = 1: length(traces(1,:))
    subplot(3,1,1)
    imagesc(max(ms.SFPs,[],3));
    title(['Cell', num2str(dup(i))])
    subplot(3,1,2)
    imagesc(footprints(:,:,i));
    subplot(3,1,3)
    plot(traces(:,i))
    pause
end

end