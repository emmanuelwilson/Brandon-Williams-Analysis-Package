function getPathROIs(paths)
    clc
    fprintf(['Define Path ROIs:\n']);
    
    cornerROIs = nan(4,2);
    for p = paths'
        s = load(p{1});
        
        if isfield(s.processed,'roi')
            continue
        end
        
        fprintf(['\t' num2str(p{1}) '\n'])
        
        figure(1)
        set(gcf,'position',[50 50 nanmax(s.processed.p').*10])
        plot(s.processed.p(1,:),s.processed.p(2,:));
        hold on
        set(gca,'xlim',[0 nanmax(s.processed.p(1,:))],'ylim',[0 nanmax(s.processed.p(2,:))])
        drawnow
        
        fprintf(['\t\tDraw Corner Rectangles.'])
        for k = 1:2
            fprintf(['\n\t\t\tCorner:\t' num2str(k)])
            cornerROIs(:,k) = getrect(); %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
        end
        s.processed.roi.corner = cornerROIs;
        for i = 1:length(cornerROIs(1,:))
            rectangle('position',cornerROIs(:,i),'linewidth',2);
        end
        drawnow
        close all
        drawnow
        
        save(p{1},'-struct','s','-v7.3');
    end
end