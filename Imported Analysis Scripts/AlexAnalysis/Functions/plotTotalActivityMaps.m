function plotTotalActivityMaps(paths)
    for p = paths'
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        m = mkTraceMaps(s.processed.p,s.processed.trace,[]);
        
        m(:,:,s.processed.splithalf.p > 0.05) = [];

        figure(1)
        set(gcf,'position',[50 50 350 150])
        imagesc(nanmean(m,3))
        colormap('jet')
        alpha(double(~isnan(m(:,:,1))))
        axis equal
        axis off
            
        slashInds = find(ismember(p{1},'/'));
        outP = ['Plots/TotalActivityMaps/' p{1}(slashInds+1:end-4)];
        saveFig(gcf,outP,'tiff')
        saveFig(gcf,outP,'pdf')
        close all
        drawnow
    end
end