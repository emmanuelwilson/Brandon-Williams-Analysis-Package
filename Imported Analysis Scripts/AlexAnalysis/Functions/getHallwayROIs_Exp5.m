function getPathROIs(paths)
    clc
    fprintf(['Define Hallway ROIs:\n']);
    
    doorROI = nan(4,2);
    roomROI = nan(4,1);
    for p = paths'
        s = load(p{1});
        
%         if isfield(s.processed,'roi')
%             continue
%         end
        
        fprintf(['\t' num2str(p{1}) '\n'])
        
        figure(1)
        set(gcf,'position',[50 50 nanmax(s.processed.p').*10])
        plot(s.processed.p(1,:),s.processed.p(2,:));
        hold on
        set(gca,'xlim',[0 nanmax(s.processed.p(1,:))],'ylim',[0 nanmax(s.processed.p(2,:))])
        drawnow
        
        fprintf(['\t\tSelect Hallway.'])
        rect = ginput(); 
        s.processed.roi.hallway = rect;
        
        
        isIn = inpolygon(s.processed.p(1,:)',s.processed.p(2,:)',rect(:,1),rect(:,2));
        close all
        figure(1)
        plot(s.processed.p(1,:),s.processed.p(2,:),'color',[0.5 0.5 0.5])
        hold on
        plot(s.processed.p(1,isIn),s.processed.p(2,isIn),'color','k',...
            'linestyle','none','marker','o','markersize',2)
        set(gcf,'position',[50 50 nanmax(s.processed.p').*10])
        set(gca,'xlim',[0 nanmax(s.processed.p(1,:))],'ylim',[0 nanmax(s.processed.p(2,:))])
        drawnow

        fprintf(['\tSelect Linear Path for collapse.'])
        rect = ginput(); 
        s.processed.roi.hallway_linear = rect;
               
%         cp = help_collapseToLine(s.processed.p(:,isIn),rect);
%         [maps samp] = mkLinMaps(cp,s.processed.trace(:,isIn));
        
%         drawnow
        pause(5)
        close all
        drawnow
        
        save(p{1},'-struct','s','-v7.3');
    end
end