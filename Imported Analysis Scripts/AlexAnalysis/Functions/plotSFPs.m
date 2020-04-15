function plotSFPs(paths)
    clc
    figure
    set(gcf,'position',[50 50 1800 900])
    for p = paths'
        s = load(p{1});
        fprintf(['Plotting Footprints:  ' p{1} '\n'])

        SFPs = (s.calcium.SFPs.^2);
        if isfield(s.processed,'exclude')
            SFPs = SFPs(:,:,s.processed.exclude.SFPs);
        end
        SFPs = SFPs./repmat(nanmax(nanmax(SFPs,[],1),[],2),[size(SFPs(:,:,1)) 1]);
        
        subplot(4,8,find(ismember(paths,p)));
        imagesc(nanmax(SFPs,[],3));
        axis off
        
%         tc = hsv(length(SFPs(1,1,:)));
%         SFPs = repmat(permute(tc,[3 4 2 1]),[size(SFPs(:,:,1,1)) 1 1]) .* ...
%             repmat(permute(SFPs,[1 2 4 3]),[1 1 3 1]);
%         imagesc(nanmax(SFPs,[],4).*1)
        
    end
    slashInds = find(ismember(p{1},'/'));
    outP = ['Plots/SpatialFootprints/' p{1}(slashInds+1:end-4)];
    saveFig(gcf,outP,'tiff')
    saveFig(gcf,outP,'pdf')
    close all
    drawnow
end